//! Alignment pipeline: chunked streaming FASTQ → rayon align → ordered write of
//! BAM / unmapped.fq.gz / align_qc.json.
//!
//! Engineering invariants:
//! - Determinism: parallel within a chunk, merged by chunk order; output record
//!   order == input order (byte-identical at any thread count).
//! - Constant memory: CHUNK pairs per block; stats and output advance per chunk.
//! - unmapped.fq.gz hard contract: every unplaced read (incl. single end of a
//!   PE pair) is written, name keeps the /1 /2 mate marker.
//! - BAM unmapped records written with FLAG 0x4 (+ mate flags).

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use rayon::prelude::*;

use crate::align::{AlignConfig, Aligner, CachedPlacements};
use crate::bam::{self, flag, BamRecord};
use crate::error::AlignError;
use crate::extend::CigarOp;
use crate::fastq::{FastqReader, FastqRecord, RecordSource};
use crate::gtf::{Junction, JunctionLib};
use crate::index::Index;
use crate::mapq::{mapq_of, ReadAlignment, PILEUP_GATE};
use crate::pair::{relate, template_length, PairRelation};
use crate::seed::Strand;
use crate::stats::{AlignStats, StatsAcc};

/// Preserve original-search limitation on an unplaced record for later rescue.
fn mark_unplaced_seed_limit(rec: &mut BamRecord, limited: bool) {
    if rec.aln.is_none() && limited {
        rec.unmapped_tags.push(esperanto_bamio::RawTag(
            *b"EK",
            esperanto_bamio::TagValue::Str("mapq_src=seed_budget;seed_search=limited".into()),
        ));
    }
}

/// Streaming block size (pairs / reads).
pub const CHUNK: usize = 100_000;

/// Pipeline output sinks (bam sink ownership is taken at run time —
/// multithreaded BGZF requires 'static).
pub struct PipelineOut<'a> {
    pub bam: Option<Box<dyn Write + Send>>,
    pub unmapped_fq: Box<dyn Write + Send>,
    /// Optional pass2-only candidate evidence, never an alignment source.
    pub diagnostics: Option<Box<dyn Write + Send>>,
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// L1 transcriptome-first engine bundle (None = genomic alignment only).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
}

/// pass1 junction discoveries: (junction, support, summed MAPQ of supporters).
pub type Discoveries = Vec<(Junction, u32, u64)>;

/// What pass 1 decided for one filtered discovery.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DiscoveryFate {
    /// No annotated junction within 3 bp: genuinely novel, merged.
    Kept,
    /// Within 3 bp of an annotated junction; pass-1 support within a 2:1 band
    /// of the annotated boundary's own observation count. Both stay merged.
    TieKept,
    /// Within 3 bp of an annotated junction and observed >2x as often:
    /// merged, and the annotated copy leaves the merged library (the sample
    /// does not express it as annotated).
    OutvotesAnnotated,
    /// Within 3 bp of an annotated junction observed >2x as often: the
    /// curated boundary wins; reported but not merged.
    OutvotedByAnnotated,
}

/// A pass-1 discovery that passed the quality gates: (junction, support, fate).
pub type FilteredDiscovery = (Junction, u32, DiscoveryFate);

/// 2-pass outputs: (stats, pass1 filtered discoveries, final library list,
/// count of annotated junctions removed by the support vote).
pub type TwoPassOut = (
    AlignStats,
    Vec<FilteredDiscovery>,
    Vec<(Junction, u32)>,
    usize,
);

fn make_aligner<'a>(
    index: &'a Index,
    cfg: AlignConfig,
    lib: &Option<Arc<JunctionLib>>,
    l1: &Option<Arc<esperanto_engine::L1Index>>,
) -> Aligner<'a> {
    let a = Aligner::new(index, cfg);
    let a = match lib {
        Some(l) => a.with_lib(l.clone()),
        None => a,
    };
    match l1 {
        Some(l1) => a.with_l1(l1.clone()),
        None => a,
    }
}

/// Alignment worker stack size. On repeat-dense reads the aligner's
/// extension/arbitration chain comes close to rayon's 2 MiB default worker
/// stack, and the crossbeam epoch collector riding on the same stack then
/// tips it over (observed as an abort at 128 threads on the full real
/// sample). 64 MiB is reserved address space only — committed lazily.
const WORKER_STACK: usize = 64 << 20;

/// Single-end pipeline (returns stats + pass1 discoveries, empty in DNA mode).
pub fn run_se(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<(AlignStats, Discoveries), AlignError> {
    let t0 = std::time::Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(WORKER_STACK)
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let mut sinks = crate::sinks::RecordSinks::new(out.bam.take(), out.index, threads)?;
    let mut unm = crate::sinks::UnmappedFq::new(&mut out.unmapped_fq);
    let mut stats = StatsAcc::new();
    let mut discoveries: Discoveries = Vec::new();

    let mut reader = FastqReader::open(r1)?;

    loop {
        let mut chunk = Vec::with_capacity(CHUNK);
        for _ in 0..CHUNK {
            match reader.next_record()? {
                Some(r) => chunk.push(r),
                None => break,
            }
        }
        if chunk.is_empty() {
            break;
        }
        let index = out.index;
        let cfg = out.config;
        let lib = out.jlib.clone();
        let l1_clone = out.l1.clone();
        let diagnostics_enabled = out.diagnostics.is_some();
        // align_read(seq) depends only on seq, not qual/name.
        // Within-chunk dedup: align each unique seq once, fan out to all reads
        // sharing it (each keeps its own name/qual). Deterministic algorithm +
        // ordered write ⇒ byte-identical.
        let (unique_seqs, uid_of_read) = dedup_se(&chunk);
        let cached: Vec<(
            Option<ReadAlignment>,
            Option<crate::diagnostics::CandidatePool>,
            bool,
        )> = {
            let miss_seqs: Vec<&[u8]> = unique_seqs.to_vec();
            pool.install(|| {
                let lib = lib.clone();
                miss_seqs
                    .par_iter()
                    .map_init(
                        move || make_aligner(index, cfg, &lib, &l1_clone),
                        |aligner, &seq| {
                            let (aln, cache) = aligner.align_read_full(seq);
                            let limited = cache.as_ref().is_some_and(|c| c.seed_search_limited());
                            let diagnostic = diagnostics_enabled
                                .then(|| diagnostic_pool(index, seq.len(), cache.as_ref()));
                            (aln, diagnostic, limited)
                        },
                    )
                    .collect()
            })
        };
        // Build BamRecord per read in input order, cloning the cached
        // alignment and pairing it with the read's own name/qual.
        let mut results: Vec<BamRecord> = Vec::with_capacity(chunk.len());
        for (i, rec) in chunk.iter().enumerate() {
            let mut br = bam::record_se(
                &rec.name,
                &rec.seq,
                &rec.qual,
                cached[uid_of_read[i]].0.as_ref(),
            );
            mark_unplaced_seed_limit(&mut br, cached[uid_of_read[i]].2);
            results.push(br);
        }
        for (i, (rec, br)) in chunk.iter().zip(&results).enumerate() {
            let (aln, diagnostic, _) = &cached[uid_of_read[i]];
            if let (Some(sink), Some(pool)) = (out.diagnostics.as_mut(), diagnostic) {
                crate::diagnostics::write_row(
                    &mut **sink,
                    rec,
                    None,
                    aln.is_some(),
                    pool.clone(),
                    None,
                )?;
            }
            stats.push_read(aln.is_some(), br.mapq);
            if let Some(x) = aln {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(br.mapq)));
                }
            }
            if aln.is_none() {
                unm.write(rec, None)?;
            }
            sinks.write(br)?;
        }
    }
    sinks.finish()?;
    unm.finish()?;
    if let Some(sink) = out.diagnostics.as_mut() {
        crate::diagnostics::flush(&mut **sink)?;
    }
    let mode = if out.config.rna { "rna-se" } else { "dna-se" };
    Ok((
        stats.finalize(mode, false, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// Single-end deduplication: group reads in a chunk by identical seq → unique seq
/// slices (by uid) + per-read → uid map.
fn dedup_se(chunk: &[FastqRecord]) -> (Vec<&[u8]>, Vec<usize>) {
    let mut seq_to_uid: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::with_capacity(chunk.len());
    let mut unique_seqs: Vec<&[u8]> = Vec::new();
    let mut uid_of_read: Vec<usize> = Vec::with_capacity(chunk.len());
    for r in chunk {
        let seq: &[u8] = r.seq.as_slice();
        let uid = *seq_to_uid.entry(seq).or_insert_with(|| {
            let u = unique_seqs.len();
            unique_seqs.push(seq);
            u
        });
        uid_of_read.push(uid);
    }
    (unique_seqs, uid_of_read)
}

/// Paired-end deduplication: group both mates of a chunk by identical seq → unique
/// seq slices (by uid, first-seen order) + per-pair → uid maps for R1/R2.
/// One unified table: a sequence appearing as both an R1 and an R2 aligns once.
#[cfg(test)]
fn dedup_pe(chunk: &[(FastqRecord, FastqRecord)]) -> (Vec<&[u8]>, Vec<usize>, Vec<usize>) {
    dedup_pe_fragments(chunk, &vec![None; chunk.len()])
}

fn dedup_pe_fragments<'a>(
    chunk: &'a [(FastqRecord, FastqRecord)],
    fragments: &[Option<crate::fragment::FragmentEvidence>],
) -> (Vec<&'a [u8]>, Vec<usize>, Vec<usize>) {
    assert_eq!(chunk.len(), fragments.len());
    let mut seq_to_uid: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::with_capacity(2 * chunk.len());
    let mut unique_seqs: Vec<&[u8]> = Vec::new();
    let mut uid1_of_pair: Vec<usize> = Vec::with_capacity(chunk.len());
    let mut uid2_of_pair: Vec<usize> = Vec::with_capacity(chunk.len());
    for ((a, b), fragment) in chunk.iter().zip(fragments) {
        let n1 = fragment.map_or(a.seq.len(), |f| f.fragment_length);
        let n2 = fragment.map_or(b.seq.len(), |f| f.fragment_length);
        uid1_of_pair.push(intern_seq(&mut seq_to_uid, &mut unique_seqs, &a.seq[..n1]));
        uid2_of_pair.push(intern_seq(&mut seq_to_uid, &mut unique_seqs, &b.seq[..n2]));
    }
    (unique_seqs, uid1_of_pair, uid2_of_pair)
}

/// Intern one sequence into the uid table; uid = first-seen insertion order.
fn intern_seq<'a>(
    seq_to_uid: &mut std::collections::HashMap<&'a [u8], usize>,
    unique_seqs: &mut Vec<&'a [u8]>,
    seq: &'a [u8],
) -> usize {
    let next_uid = unique_seqs.len();
    *seq_to_uid.entry(seq).or_insert_with(|| {
        unique_seqs.push(seq);
        next_uid
    })
}

fn diagnostic_pool(
    index: &Index,
    read_len: usize,
    cache: Option<&CachedPlacements>,
) -> crate::diagnostics::CandidatePool {
    match cache {
        Some(cache) => cache.diagnostic_pool(index, read_len, "independent"),
        None => crate::diagnostics::CandidatePool::capture(index, read_len, &[], "independent"),
    }
}

/// Paired-end pipeline (returns stats + pass1 discoveries, empty in DNA mode).
pub fn run_pe(
    out: &mut PipelineOut,
    r1: &Path,
    r2: &Path,
    threads: usize,
) -> Result<(AlignStats, Discoveries), AlignError> {
    let t0 = std::time::Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(WORKER_STACK)
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let mut sinks = crate::sinks::RecordSinks::new(out.bam.take(), out.index, threads)?;
    let mut unm = crate::sinks::UnmappedFq::new(&mut out.unmapped_fq);
    let mut stats = StatsAcc::new();
    let mut discoveries: Discoveries = Vec::new();
    let mut rd1 = FastqReader::open(r1)?;
    let mut rd2 = FastqReader::open(r2)?;

    loop {
        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK);
        for _ in 0..CHUNK {
            match (rd1.next_record()?, rd2.next_record()?) {
                (Some(a), Some(b)) => chunk.push((a, b)),
                (None, None) => break,
                // Pair sync: either end hitting EOF first is an error.
                (a, b) => {
                    return Err(AlignError::FastqFormat {
                        line: 0,
                        msg: format!(
                            "R1/R2 record count mismatch (r1 has more: {}, r2 has more: {})",
                            a.is_some(),
                            b.is_some()
                        ),
                    });
                }
            }
        }
        if chunk.is_empty() {
            break;
        }
        let index = out.index;
        let cfg = out.config;
        let lib = out.jlib.clone();
        let l1_clone = out.l1.clone();
        // Paired-end deduplication — two layers, chunk-local only (the table dies
        // with the chunk to keep memory bounded; no cross-chunk cache):
        //   1. primary alignment per UNIQUE effective sequence. Raw pair
        //      endpoint evidence is computed before interning; exact effective
        //      bytes are the key, so different mate support cannot contaminate
        //      a sequence-only cache. Full effective seq bytes are the
        //      key: name/quality never participate in primary alignment
        //      today — if alignment ever consumes quality, this dedup must
        //      be revisited (the key must then include quality).
        //   2. per-pair mate-rescue arbitration is NEVER cached (rescue
        //      depends on the concrete pair): each pair clones the cached
        //      primary of its two sequences and runs the existing rescue
        //      logic unchanged.
        let fragments: Vec<_> = chunk
            .iter()
            .map(|(a, b)| {
                (cfg.rna && cfg.fragment_endpoints)
                    .then(|| crate::fragment::infer(&a.seq, &b.seq))
                    .flatten()
            })
            .collect();
        let (unique_seqs, uid1_of_pair, uid2_of_pair) = dedup_pe_fragments(&chunk, &fragments);
        let unique_alns: Vec<(Option<ReadAlignment>, Option<CachedPlacements>)> =
            pool.install(|| {
                let lib = lib.clone();
                let l1_clone = l1_clone.clone();
                unique_seqs
                    .par_iter()
                    .map_init(
                        move || make_aligner(index, cfg, &lib, &l1_clone),
                        |aligner, &seq| aligner.align_read_full(seq),
                    )
                    .collect()
            });
        // Parallel per-pair rescue arbitration reads the uid-indexed primary
        // results (order preserved by collect); BAM records are built in the
        // serial write section.
        let diagnostics_enabled = out.diagnostics.is_some();
        type PairResult = (
            Option<ReadAlignment>,
            Option<ReadAlignment>,
            Option<crate::diagnostics::CandidatePool>,
            Option<crate::diagnostics::CandidatePool>,
        );
        let results: Vec<PairResult> = pool.install(|| {
            let lib = lib.clone();
            let l1_clone = l1_clone.clone();
            chunk
                .par_iter()
                .enumerate()
                .map_init(
                    move || make_aligner(index, cfg, &lib, &l1_clone),
                    |aligner, (i, (a, b))| {
                        let seq1 = unique_seqs[uid1_of_pair[i]];
                        let seq2 = unique_seqs[uid2_of_pair[i]];
                        let (base1, cache1) = &unique_alns[uid1_of_pair[i]];
                        let (base2, cache2) = &unique_alns[uid2_of_pair[i]];
                        let mut r1 = base1.clone();
                        let mut r2 = base2.clone();
                        let mut diag1 = None;
                        let mut diag2 = None;
                        // Pair-anchored arbitration: a mate whose own
                        // placement is poor (unplaced, mostly clipped, or
                        // below the pileup gate) is re-arbitrated with the
                        // locus its confidently placed mate points at as an
                        // extra candidate (see `Aligner::align_read_with_mate`).
                        let confident = |r: &Option<ReadAlignment>| {
                            r.as_ref()
                                .is_some_and(|m| !m.rescued && mapq_of(m) >= PILEUP_GATE)
                        };
                        if placement_poor(&r1, seq1.len()) && confident(&r2) {
                            let (mc, mp) = r2.as_ref().map(|m| (m.contig, m.pos)).unwrap();
                            if let Some(n) = aligner.align_read_with_mate_observed(
                                seq1,
                                mc,
                                mp,
                                cache1.as_ref(),
                                diagnostics_enabled.then_some(&mut diag1),
                            ) {
                                r1 = Some(n);
                            }
                        }
                        if placement_poor(&r2, seq2.len()) && confident(&r1) {
                            let (mc, mp) = r1.as_ref().map(|m| (m.contig, m.pos)).unwrap();
                            if let Some(n) = aligner.align_read_with_mate_observed(
                                seq2,
                                mc,
                                mp,
                                cache2.as_ref(),
                                diagnostics_enabled.then_some(&mut diag2),
                            ) {
                                r2 = Some(n);
                            }
                        }
                        // Qualification is pair-specific and never interned by
                        // sequence. Already placed records keep their old result.
                        if r1.is_none() && confident(&r2) {
                            r1 = aligner.recover_short_body_observed(
                                seq1,
                                r2.as_ref().unwrap(),
                                cache1.as_ref(),
                                diagnostics_enabled.then_some(&mut diag1),
                            );
                        }
                        if r2.is_none() && confident(&r1) {
                            r2 = aligner.recover_short_body_observed(
                                seq2,
                                r1.as_ref().unwrap(),
                                cache2.as_ref(),
                                diagnostics_enabled.then_some(&mut diag2),
                            );
                        }
                        if cfg.rna && cfg.pair_substitution && fragments[i].is_some() {
                            if let (Some(old1), Some(old2), Some(c1), Some(c2)) =
                                (&r1, &r2, cache1.as_ref(), cache2.as_ref())
                            {
                                if let Some(witness) = crate::align::pair_substitution::witness(
                                    index,
                                    [seq1, seq2],
                                    [&a.qual[..seq1.len()], &b.qual[..seq2.len()]],
                                    [old1, old2],
                                ) {
                                    let next1 = aligner.extend_pair_substitution_observed(
                                        seq1,
                                        old1,
                                        c1,
                                        &witness,
                                        diagnostics_enabled.then_some(&mut diag1),
                                    );
                                    let next2 = aligner.extend_pair_substitution_observed(
                                        seq2,
                                        old2,
                                        c2,
                                        &witness,
                                        diagnostics_enabled.then_some(&mut diag2),
                                    );
                                    // Both decisions used the unchanged original pair.
                                    if let Some(next) = next1 {
                                        r1 = Some(next);
                                    }
                                    if let Some(next) = next2 {
                                        r2 = Some(next);
                                    }
                                }
                            }
                        }
                        for (aln, raw_len, effective_len) in [
                            (&mut r1, a.seq.len(), seq1.len()),
                            (&mut r2, b.seq.len(), seq2.len()),
                        ] {
                            if let Some(aln) = aln {
                                crate::fragment::lift_cigar(
                                    &mut aln.cigar,
                                    aln.strand,
                                    raw_len - effective_len,
                                );
                            }
                        }
                        (r1, r2, diag1, diag2)
                    },
                )
                .collect()
        });
        for (i, ((a, b), (aln1, aln2, diag1, diag2))) in chunk.iter().zip(&results).enumerate() {
            if let Some(sink) = out.diagnostics.as_mut() {
                for (read, mate, aln, uid, attempt) in [
                    (a, 1, aln1, uid1_of_pair[i], diag1),
                    (b, 2, aln2, uid2_of_pair[i], diag2),
                ] {
                    let effective_len = unique_seqs[uid].len();
                    crate::diagnostics::write_row_with_fragment(
                        &mut **sink,
                        read,
                        Some(mate),
                        aln.is_some(),
                        diagnostic_pool(index, effective_len, unique_alns[uid].1.as_ref()),
                        attempt.clone(),
                        fragments[i],
                    )?;
                }
            }
            // Records carry every number the stats/discovery section needs
            // (MAPQ, proper-pair flag, tlen) — one derivation per field,
            // no re-computation from the alignment.
            let (mut ra, mut rb) = pair_records(
                &a.name, &a.seq, &a.qual, aln1, &b.name, &b.seq, &b.qual, aln2,
            );
            if let Some(fragment) = fragments[i] {
                for record in [&mut ra, &mut rb] {
                    let tags = match &mut record.aln {
                        Some(aln) => &mut aln.tags,
                        None => &mut record.unmapped_tags,
                    };
                    tags.push(esperanto_bamio::RawTag(
                        *b"XF",
                        esperanto_bamio::TagValue::Int(fragment.fragment_length as i32),
                    ));
                }
            }
            mark_unplaced_seed_limit(
                &mut ra,
                unique_alns[uid1_of_pair[i]]
                    .1
                    .as_ref()
                    .is_some_and(|c| c.seed_search_limited()),
            );
            mark_unplaced_seed_limit(
                &mut rb,
                unique_alns[uid2_of_pair[i]]
                    .1
                    .as_ref()
                    .is_some_and(|c| c.seed_search_limited()),
            );
            if let Some(x) = aln1 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(ra.mapq)));
                }
            }
            if let Some(x) = aln2 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(rb.mapq)));
                }
            }
            stats.push_read(aln1.is_some(), ra.mapq);
            stats.push_read(aln2.is_some(), rb.mapq);
            if let Some(x) = aln1 {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if let Some(x) = aln2 {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if aln1.is_none() {
                if aln2.is_some() {
                    stats.push_rescue_fail();
                }
                unm.write(a, Some(1))?;
            }
            if aln2.is_none() {
                if aln1.is_some() {
                    stats.push_rescue_fail();
                }
                unm.write(b, Some(2))?;
            }
            // proper ⇔ both placed and relate()==Proper (mk set the flag from
            // exactly that condition); tlen is the value mk stored for R1.
            if ra.flag & flag::PROPER_PAIR != 0 {
                stats.push_proper_pair();
                if let Some((_, _, tlen)) = ra.mate {
                    stats.push_insert(tlen);
                }
            }
            sinks.write(&ra)?;
            sinks.write(&rb)?;
        }
    }
    sinks.finish()?;
    unm.finish()?;
    if let Some(sink) = out.diagnostics.as_mut() {
        crate::diagnostics::flush(&mut **sink)?;
    }
    let mode = if out.config.rna { "rna-pe" } else { "dna-pe" };
    Ok((
        stats.finalize(mode, true, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// Build PE paired records (FLAG/mate fields per spec).
///
/// Alignments are borrowed: no per-pair `ReadAlignment` clone — the records
/// copy out every field they need.
#[allow(clippy::too_many_arguments)]
fn pair_records(
    n1: &[u8],
    s1: &[u8],
    q1: &[u8],
    a1: &Option<ReadAlignment>,
    n2: &[u8],
    s2: &[u8],
    q2: &[u8],
    a2: &Option<ReadAlignment>,
) -> (BamRecord, BamRecord) {
    let proper = match (a1, a2) {
        (Some(x), Some(y)) => relate(x, y) == PairRelation::Proper,
        _ => false,
    };
    let tlen = match (a1, a2) {
        (Some(x), Some(y)) if proper => template_length(x, y),
        _ => 0,
    };
    let mk = |name: &[u8],
              seq: &[u8],
              qual: &[u8],
              aln: &Option<ReadAlignment>,
              mate_aln: &Option<ReadAlignment>,
              self_bit: u16,
              tlen: i32| {
        let mut f = flag::PAIRED | self_bit;
        if proper {
            f |= flag::PROPER_PAIR;
        }
        let (mate_c, mate_p) = match mate_aln {
            Some(m) => {
                if m.strand == Strand::Minus {
                    f |= flag::MATE_REVERSE;
                }
                (m.contig as i32, m.pos as i32)
            }
            None => {
                f |= flag::MATE_UNMAPPED;
                (-1, -1)
            }
        };
        let q = match aln {
            Some(a) => {
                if a.strand == Strand::Minus {
                    f |= flag::REVERSE;
                }
                mapq_of(a)
            }
            None => {
                f |= flag::UNMAPPED;
                0
            }
        };
        // Reference-forward SEQ direction (same as record_se): minus strand stores the
        // reference-forward SEQ (revcomp) with reversed QUAL; unmapped keeps
        // the original orientation.
        let (s_out, q_out) = crate::bam::t13_phred(
            aln.as_ref().is_some_and(|a| a.strand == Strand::Minus),
            seq,
            qual,
        );
        BamRecord {
            name: name.to_vec(),
            flag: f,
            aln: aln.as_ref().map(crate::bam::aln_view),
            unmapped_tags: Vec::new(),
            mapq: q,
            seq: s_out,
            qual: q_out,
            mate: Some((mate_c, mate_p, tlen)),
        }
    };
    (
        mk(n1, s1, q1, a1, a2, flag::READ1, tlen),
        mk(n2, s2, q2, a2, a1, flag::READ2, -tlen),
    )
}

/// SE 2-pass (same semantics as the PE version). Runs pass-1 discovery even
/// with the L1 engine on: the annotation covers only ~75% of the junctions
/// expressed in a typical sample, and the L1 path cannot place a read across
/// a junction neither the bundle nor the library knows.
pub fn run_se_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    let base_lib = out.jlib.take();
    out.jlib = base_lib.clone();
    let real_bam = out.bam.take();
    out.bam = Some(Box::new(std::io::sink()));
    let real_unm = std::mem::replace(&mut out.unmapped_fq, Box::new(std::io::sink()));
    let real_diagnostics = out.diagnostics.take();
    let pass1 = run_se(out, r1, threads);
    let sink_unm = std::mem::replace(&mut out.unmapped_fq, real_unm);
    drop(sink_unm);
    out.bam = real_bam;
    out.diagnostics = real_diagnostics;
    let (_, discovered) = pass1?;
    // Merge: full sjdb + discoveries with ≥2 support.
    let mut merged: Vec<(Junction, u32)> = match &base_lib {
        Some(lib) => lib
            .junctions
            .iter()
            .copied()
            .zip(lib.counts.iter().copied())
            .collect(),
        None => Vec::new(),
    };
    let (kept, removed_annotated) = filter_discoveries(discovered, base_lib.as_deref());
    if !removed_annotated.is_empty() {
        let removed: std::collections::BTreeSet<&Junction> = removed_annotated.iter().collect();
        merged.retain(|(j, _)| !removed.contains(j));
    }
    merged.extend(
        kept.iter()
            .filter(|(_, _, fate)| *fate != DiscoveryFate::OutvotedByAnnotated)
            .map(|(j, c, _)| (*j, *c)),
    );
    let final_lib = JunctionLib::build_with_counts(merged);
    let final_list: Vec<(Junction, u32)> = final_lib
        .junctions
        .iter()
        .copied()
        .zip(final_lib.counts.iter().copied())
        .collect();
    out.jlib = Some(Arc::new(final_lib));
    let (stats, _) = run_se(out, r1, threads)?;
    Ok((stats, kept, final_list, removed_annotated.len()))
}

/// RNA 2-pass: pass1 discovers junctions (>=2 support into the library) →
/// merge with sjdb → pass2 aligns with the merged library. Runs even with
/// the L1 engine on: the annotation covers only ~75% of the junctions
/// expressed in a typical sample, and the L1 path cannot place a read across
/// a junction neither the bundle nor the library knows.
pub fn run_pe_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    r2: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    // pass1: discovery round (BAM/unmapped discarded; real artifacts only
    // come from pass2).
    let base_lib = out.jlib.take();
    out.jlib = base_lib.clone();
    let real_bam = out.bam.take();
    out.bam = Some(Box::new(std::io::sink()));
    let real_unm = std::mem::replace(&mut out.unmapped_fq, Box::new(std::io::sink()));
    let real_diagnostics = out.diagnostics.take();
    // Short-body MAPQ0 qualification cannot change the discovery library and
    // thereby perturb ordinary pass2 placements. It runs only in final output.
    let mate_short_body = out.config.mate_short_body;
    let pair_substitution = out.config.pair_substitution;
    out.config.mate_short_body = false;
    out.config.pair_substitution = false;
    let pass1 = run_pe(out, r1, r2, threads);
    out.config.mate_short_body = mate_short_body;
    out.config.pair_substitution = pair_substitution;
    let sink_unm = std::mem::replace(&mut out.unmapped_fq, real_unm);
    drop(sink_unm);
    out.bam = real_bam;
    out.diagnostics = real_diagnostics;
    let (_, discovered) = pass1?;
    let mut merged: Vec<(Junction, u32)> = match &base_lib {
        Some(lib) => lib
            .junctions
            .iter()
            .copied()
            .zip(lib.counts.iter().copied())
            .collect(),
        None => Vec::new(),
    };
    let (kept, removed_annotated) = filter_discoveries(discovered, base_lib.as_deref());
    if !removed_annotated.is_empty() {
        let removed: std::collections::BTreeSet<&Junction> = removed_annotated.iter().collect();
        merged.retain(|(j, _)| !removed.contains(j));
    }
    merged.extend(
        kept.iter()
            .filter(|(_, _, fate)| *fate != DiscoveryFate::OutvotedByAnnotated)
            .map(|(j, c, _)| (*j, *c)),
    );
    let final_lib = JunctionLib::build_with_counts(merged);
    let final_list: Vec<(Junction, u32)> = final_lib
        .junctions
        .iter()
        .copied()
        .zip(final_lib.counts.iter().copied())
        .collect();
    out.jlib = Some(Arc::new(final_lib));
    // pass2: production alignment with the merged library.
    let (stats, _) = run_pe(out, r1, r2, threads)?;
    Ok((stats, kept, final_list, removed_annotated.len()))
}

/// Junction-level filter for 2-pass: support >= 2 AND mean supporter MAPQ >= 20
/// AND span <= 500kb (low-quality cross-repeat-family junctions stay out).
///
/// A surviving discovery within 3 bp of an annotated junction (same contig,
/// either strand) is a shifted re-observation of a curated boundary. Keeping
/// both creates a same-score tie that flips boundaries at scale; always
/// dropping the discovery loses real signal where the annotation is the side
/// the sample does not express. So the pass-1 observation counts vote: the
/// side observed >2x as often stays and the other leaves the merged library;
/// inside a 1:2..2:1 band both stay and arbitration decides. The annotated
/// side's pass-1 support is its own exact re-observation count in the
/// discovery stream (0 when never observed). The strand check is deliberately
/// blind: a shifted copy can surface under the opposite strand label
/// (observed on sim1).
pub fn filter_discoveries(
    discovered: Discoveries,
    annotated: Option<&JunctionLib>,
) -> (Vec<FilteredDiscovery>, Vec<Junction>) {
    let mut acc: std::collections::BTreeMap<Junction, (u32, u64)> =
        std::collections::BTreeMap::new();
    for (j, c, mq) in discovered {
        let e = acc.entry(j).or_default();
        e.0 += c;
        e.1 += mq;
    }
    let mut out = Vec::new();
    let mut removed: std::collections::BTreeSet<Junction> = std::collections::BTreeSet::new();
    for (j, (c, mqsum)) in &acc {
        if !(*c >= 2 && *mqsum >= u64::from(*c) * 20 && j.end.saturating_sub(j.start) <= 500_000) {
            continue;
        }
        let fate = match annotated {
            None => DiscoveryFate::Kept,
            Some(lib) => {
                let (js, _) = lib.range_start(
                    j.contig,
                    j.start.saturating_sub(3),
                    j.start.saturating_add(4),
                );
                let nbrs: Vec<&Junction> = js
                    .iter()
                    .filter(|a| a.start.abs_diff(j.start) <= 3 && a.end.abs_diff(j.end) <= 3)
                    .collect();
                if nbrs.is_empty() {
                    DiscoveryFate::Kept
                } else {
                    let sup_a = nbrs
                        .iter()
                        .map(|a| acc.get(a).map_or(0, |(s, _)| *s))
                        .max()
                        .unwrap_or(0);
                    if *c > 2 * sup_a {
                        for a in nbrs {
                            if *c > 2 * acc.get(a).map_or(0, |(s, _)| *s) {
                                removed.insert(*a);
                            }
                        }
                        DiscoveryFate::OutvotesAnnotated
                    } else if sup_a > 2 * *c {
                        DiscoveryFate::OutvotedByAnnotated
                    } else {
                        DiscoveryFate::TieKept
                    }
                }
            }
        };
        out.push((*j, *c, fate));
    }
    (out, removed.into_iter().collect())
}

/// A placement the pair-anchored arbitration should second-guess: absent,
/// explaining less than 90% of the read, or below the pileup gate.
fn placement_poor(r: &Option<ReadAlignment>, read_len: usize) -> bool {
    match r {
        None => true,
        Some(a) => {
            let aligned: u32 = a
                .cigar
                .iter()
                .map(|op| match op {
                    CigarOp::Match(n) | CigarOp::Ins(n) => *n,
                    _ => 0,
                })
                .sum();
            (aligned as f64) < 0.9 * read_len as f64 || mapq_of(a) < PILEUP_GATE
        }
    }
}

#[cfg(test)]
mod fragment_integration;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paired_fragment_pipeline_preserves_raw_records_and_recovers_short_body() {
        let body = crate::fragment::tests::body(80);
        let [a, b] = crate::fragment::tests::pair(&body, [150, 150], 0);
        let mut genome = vec![b'N'; 2080];
        genome[1000..1080].copy_from_slice(&body);
        let reference =
            crate::fasta::parse_fasta_bytes(&[b">chrT\n".as_slice(), &genome, b"\n"].concat())
                .unwrap();
        let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
        let dir =
            std::env::temp_dir().join(format!("esp-paired-fragment-red-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut inputs = Vec::new();
        let originals = [a, b];
        let raw_quals: Vec<Vec<u8>> = (0..2)
            .map(|mate| {
                (0..150)
                    .map(|i| b'!' + ((i + 7 * mate) % 41) as u8)
                    .collect()
            })
            .collect();
        for (i, seq) in originals.iter().enumerate() {
            let path = dir.join(format!("r{}.fq", i + 1));
            std::fs::write(
                &path,
                [
                    format!("@fragment/{}\n", i + 1).as_bytes(),
                    seq,
                    b"\n+\n",
                    &raw_quals[i],
                    b"\n",
                ]
                .concat(),
            )
            .unwrap();
            inputs.push(path);
        }
        let run = |paths: &[std::path::PathBuf],
                   enabled: bool,
                   rna: bool,
                   diagnostics: bool,
                   threads: usize| {
            let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
            let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut config = AlignConfig::rna_default();
            config.fragment_endpoints = enabled;
            config.rna = rna;
            let mut out = PipelineOut {
                bam: Some(Box::new(SharedSink(raw.clone()))),
                unmapped_fq: Box::new(std::io::sink()),
                diagnostics: diagnostics
                    .then(|| Box::new(SharedSink(diag.clone())) as Box<dyn Write + Send>),
                index: &index,
                config,
                jlib: None,
                l1: None,
            };
            let (stats, _) = run_pe(&mut out, &paths[0], &paths[1], threads).unwrap();
            let raw = raw.lock().unwrap().clone();
            let diag = diag.lock().unwrap().clone();
            (stats.mapped_reads, raw, diag)
        };
        let on = run(&inputs, true, true, true, 1);
        let off = run(&inputs, false, true, true, 1);
        let dna = run(&inputs, true, false, true, 1);
        let silent = run(&inputs, true, true, false, 1);
        let parallel = run(&inputs, true, true, true, 2);
        assert_eq!(
            on.0, 2,
            "raw candidates need supported fragment eligibility"
        );
        assert_eq!(off.0, 0);
        assert_eq!(dna.0, 0);
        assert_eq!(on.1, silent.1, "diagnostics cannot change primary bytes");
        assert_eq!(on.2, parallel.2, "thread-independent diagnostic rows");
        let rows: Vec<serde_json::Value> = std::str::from_utf8(&on.2)
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        assert_eq!(rows.len(), 2);
        for (i, row) in rows.into_iter().enumerate() {
            assert_eq!(row["schema_version"], 2);
            assert_eq!(row["read_length"], 150);
            assert_eq!(row["effective_read_length"], 80);
            assert_eq!(row["minimum_primary_query_bases"], 99);
            assert_eq!(row["effective_minimum_primary_query_bases"], 53);
            let winner = &row["independent"]["candidates"][0];
            assert_eq!(winner["common_score"], 160);
            assert_eq!(winner["informative_query_bases"], 80);
            assert_eq!(winner["cigar"], if i == 0 { "80M70S" } else { "70S80M" });
        }
        let mut reader = noodles::bam::io::Reader::new(std::io::Cursor::new(&on.1));
        let _ = reader.read_header().unwrap();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), 2);
        for (i, record) in records.iter().enumerate() {
            assert!(!record.flags().is_unmapped());
            assert_eq!(record.sequence().len(), 150);
            let stored: Vec<_> = record.sequence().iter().collect();
            let forward = if record.flags().is_reverse_complemented() {
                crate::align::revcomp(&stored)
            } else {
                stored
            };
            assert_eq!(forward, originals[i]);
            let mut quality: Vec<_> = record.quality_scores().iter().collect();
            if record.flags().is_reverse_complemented() {
                quality.reverse();
            }
            assert_eq!(
                quality,
                raw_quals[i].iter().map(|q| q - 33).collect::<Vec<_>>()
            );
            assert!(format!("{:?}", record.data()).contains("XF"));
        }
        // The same durable fragment boundary must survive an unplaced result.
        // Otherwise later rescue/restart would align the raw adapters again.
        let missing =
            crate::fragment::tests::pair(&crate::fragment::tests::body(300)[220..], [150, 150], 0);
        let mut missing_paths = Vec::new();
        for (i, seq) in missing.iter().enumerate() {
            let path = dir.join(format!("missing-{i}.fq"));
            std::fs::write(
                &path,
                [
                    b"@missing\n".as_slice(),
                    seq,
                    b"\n+\n",
                    &vec![b'I'; 150],
                    b"\n",
                ]
                .concat(),
            )
            .unwrap();
            missing_paths.push(path);
        }
        let absent = run(&missing_paths, true, true, true, 1);
        assert_eq!(absent.0, 0);
        let mut reader = noodles::bam::io::Reader::new(std::io::Cursor::new(&absent.1));
        reader.read_header().unwrap();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), 2);
        for record in records {
            assert!(record.flags().is_unmapped());
            assert!(
                format!("{:?}", record.data()).contains("XF"),
                "unmapped fragment boundary must be durable"
            );
        }
    }

    #[test]
    fn paired_fragment_cache_uses_effective_bytes_and_never_raw_pair_aliases() {
        let body = crate::fragment::tests::body(80);
        let [a, b] = crate::fragment::tests::pair(&body, [150, 150], 0);
        let rec = |name: &str, seq: Vec<u8>| FastqRecord {
            name: name.as_bytes().to_vec(),
            qual: vec![b'I'; seq.len()],
            seq,
        };
        let chunk = vec![
            (rec("supported/1", a.clone()), rec("supported/2", b.clone())),
            (
                rec("unsupported/1", a.clone()),
                rec("unsupported/2", vec![b'G'; 150]),
            ),
            (rec("duplicate/1", a), rec("duplicate/2", b)),
            {
                let [a, b] = crate::fragment::tests::pair(&body, [170, 150], 0);
                (rec("long/1", a), rec("long/2", b))
            },
            (
                rec("plain/1", body.clone()),
                rec("plain/2", crate::align::revcomp(&body)),
            ),
        ];
        let evidence: Vec<_> = chunk
            .iter()
            .map(|(a, b)| crate::fragment::infer(&a.seq, &b.seq))
            .collect();
        assert!(evidence[0].is_some() && evidence[1].is_none() && evidence[2].is_some());
        let (unique, uid1, uid2) = dedup_pe_fragments(&chunk, &evidence);
        assert_eq!(unique[uid1[0]], &body);
        assert_eq!(unique[uid1[1]].len(), 150);
        assert_ne!(uid1[0], uid1[1]);
        assert_eq!(uid1[0], uid1[2]);
        assert_eq!(uid2[0], uid2[2]);
        assert_eq!(
            uid1[0], uid1[3],
            "different raw lengths share effective bytes"
        );
        assert_eq!(
            uid1[0], uid1[4],
            "untrimmed read may share exact prefix bytes"
        );
        assert_eq!(uid2[0], uid2[3]);
        assert_eq!(uid2[0], uid2[4]);
    }

    #[test]
    fn paired_fragment_biological_adapter_mimic_exposes_real_body_loss() {
        let [a, b] = crate::fragment::tests::pair(&crate::fragment::tests::body(80), [150, 150], 0);
        // Alternative biological truth: a long molecule has these complete
        // reads at its two ends; the adapter-like sequence is genuine body.
        let molecule = [
            a.clone(),
            crate::fragment::tests::body(97),
            crate::align::revcomp(&b),
        ]
        .concat();
        let reference = crate::fasta::parse_fasta_bytes(
            &[b">biological\n".as_slice(), &molecule, b"\n"].concat(),
        )
        .unwrap();
        let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        let evidence = crate::fragment::infer(&a, &b).unwrap();
        assert_eq!(evidence.fragment_length, 80);
        let mut clipped_real_bases = 0;
        for raw in [&a, &b] {
            let (full, _) = aligner.align_read_full(raw);
            let full = full.expect("the complete biological read exists in reference");
            assert_eq!(full.cigar, vec![crate::extend::CigarOp::Match(150)]);
            let (short, _) = aligner.align_read_full(&raw[..80]);
            let mut short = short.unwrap();
            crate::fragment::lift_cigar(&mut short.cigar, short.strand, 70);
            clipped_real_bases += short
                .cigar
                .iter()
                .map(|op| match op {
                    crate::extend::CigarOp::SoftClip(n) => *n,
                    _ => 0,
                })
                .sum::<u32>();
        }
        // This is a real false-trimming limit, NOT an expected abstention or
        // successful recovery. Identical FASTQ bytes do not identify origin.
        assert_eq!(clipped_real_bases, 140);
    }

    #[test]
    fn paired_fragment_effective_repeat_budget_remains_mapq_zero() {
        let body = crate::fragment::tests::body(80);
        let genome = [vec![b'N'; 50], body.clone(), vec![b'N'; 50]]
            .concat()
            .repeat(2000);
        let reference =
            crate::fasta::parse_fasta_bytes(&[b">repeat\n".as_slice(), &genome, b"\n"].concat())
                .unwrap();
        let mut index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
        index.freq_cutoff = 1;
        let dir = std::env::temp_dir().join(format!("esp-fragment-budget-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let raw_pair = crate::fragment::tests::pair(&body, [150, 150], 0);
        let mut paths = Vec::new();
        for (i, seq) in raw_pair.iter().enumerate() {
            let path = dir.join(format!("r{i}.fq"));
            std::fs::write(
                &path,
                [
                    b"@repeat\n".as_slice(),
                    seq,
                    b"\n+\n",
                    &vec![b'I'; 150],
                    b"\n",
                ]
                .concat(),
            )
            .unwrap();
            paths.push(path);
        }
        let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
        let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut out = PipelineOut {
            bam: Some(Box::new(SharedSink(raw.clone()))),
            unmapped_fq: Box::new(std::io::sink()),
            diagnostics: Some(Box::new(SharedSink(diag.clone()))),
            index: &index,
            config: AlignConfig::rna_default(),
            jlib: None,
            l1: None,
        };
        let (stats, discoveries) = run_pe(&mut out, &paths[0], &paths[1], 2).unwrap();
        assert_eq!(stats.mapped_reads, 2);
        assert!(
            discoveries.is_empty(),
            "MAPQ-zero repeats cannot seed junction confidence"
        );
        let bytes = raw.lock().unwrap().clone();
        let mut reader = noodles::bam::io::Reader::new(std::io::Cursor::new(bytes));
        reader.read_header().unwrap();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), 2);
        for record in records {
            assert!(!record.flags().is_unmapped());
            assert_eq!(record.mapping_quality().unwrap().get(), 0);
            let tags = format!("{:?}", record.data());
            assert!(tags.contains("XF"));
            assert!(tags.contains("seed_search=limited"), "{tags}");
        }
        let bytes = diag.lock().unwrap();
        for line in std::str::from_utf8(&bytes).unwrap().lines() {
            let row: serde_json::Value = serde_json::from_str(line).unwrap();
            assert_eq!(row["schema_version"], 2);
            let budget = &row["independent"]["seed_budget"];
            assert!(budget["eligible_anchors"].as_u64().unwrap() > 20000);
            assert!(budget["retained_anchors"].as_u64().unwrap() <= 20000);
        }
    }

    #[test]
    fn candidate_diagnostics_are_pass2_only_ordered_and_do_not_change_primary_outputs() {
        let mut state = 0xc7735ab3u64;
        let seq: Vec<u8> = (0..4000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let reference =
            crate::fasta::parse_fasta_bytes(&[b">chrT\n".as_slice(), &seq, b"\n"].concat())
                .unwrap();
        let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
        let dir =
            std::env::temp_dir().join(format!("esp-candidate-diagnostics-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (r1, r2) = (dir.join("r1.fq"), dir.join("r2.fq"));
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (i, covered) in [96, 98, 99, 98].into_iter().enumerate() {
            // Duplicate 98-base sequences must retain distinct read/mate identities.
            let mut read = seq[1000..1000 + covered].to_vec();
            read.extend(vec![b'N'; 150 - covered]);
            let mut reverse = read.clone();
            reverse.reverse();
            for x in &mut reverse {
                *x = match *x {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                };
            }
            for (sink, name, query) in [
                (&mut a, format!("pair{i}/1"), &read),
                (&mut b, format!("pair{i}/2"), &reverse),
            ] {
                sink.extend_from_slice(format!("@{name}\n").as_bytes());
                sink.extend_from_slice(query);
                sink.extend_from_slice(b"\n+\n");
                sink.extend(vec![b'I'; 150]);
                sink.push(b'\n');
            }
        }
        std::fs::write(&r1, a).unwrap();
        std::fs::write(&r2, b).unwrap();
        let run = |enabled: bool, threads: usize| {
            let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
            let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
            let unm = Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut out = PipelineOut {
                bam: Some(Box::new(SharedSink(raw.clone()))),
                unmapped_fq: Box::new(SharedSink(unm.clone())),
                diagnostics: enabled
                    .then(|| Box::new(SharedSink(diag.clone())) as Box<dyn Write + Send>),
                index: &index,
                config: AlignConfig::rna_default(),
                jlib: None,
                l1: None,
            };
            let (stats, discoveries, library, removed) =
                run_pe_2pass(&mut out, &r1, &r2, threads).unwrap();
            drop(out);
            let result = (
                serde_json::to_value(stats).unwrap(),
                discoveries,
                library,
                removed,
                raw.lock().unwrap().clone(),
                unm.lock().unwrap().clone(),
                diag.lock().unwrap().clone(),
            );
            result
        };
        let off = run(false, 1);
        let on = run(true, 1);
        let threaded = run(true, 2);
        assert_eq!(off.0["mapped_reads"], on.0["mapped_reads"]);
        assert_eq!(off.1, on.1);
        assert_eq!(off.2, on.2);
        assert_eq!(
            off.4, on.4,
            "diagnostics must not change BAM bytes at fixed thread count"
        );
        assert_eq!(off.5, on.5, "diagnostics must not change unmapped FASTQ");
        assert_eq!(
            on.6, threaded.6,
            "JSONL must be deterministic across threads"
        );
        let rows: Vec<serde_json::Value> = std::str::from_utf8(&on.6)
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        assert_eq!(rows.len(), 8, "discovery pass must not duplicate rows");
        for (i, row) in rows.iter().enumerate() {
            let covered = [96, 98, 99, 98][i / 2];
            assert_eq!(row["read_name"], format!("pair{}/{}", i / 2, i % 2 + 1));
            assert_eq!(row["mate"], i % 2 + 1);
            assert_eq!(row["pass2_primary_placed"], covered == 99);
            assert_eq!(row["minimum_primary_query_bases"], 99);
            assert!(row["independent"]["candidates"]
                .as_array()
                .unwrap()
                .iter()
                .any(|c| c["informative_query_bases"] == covered));
        }
        let run_single = |enabled: bool| {
            let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
            let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
            let unm = Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut out = PipelineOut {
                bam: Some(Box::new(SharedSink(raw.clone()))),
                unmapped_fq: Box::new(SharedSink(unm.clone())),
                diagnostics: enabled
                    .then(|| Box::new(SharedSink(diag.clone())) as Box<dyn Write + Send>),
                index: &index,
                config: AlignConfig::rna_default(),
                jlib: None,
                l1: None,
            };
            let (_, discoveries, _, _) = run_se_2pass(&mut out, &r1, 2).unwrap();
            drop(out);
            let result = (
                discoveries,
                raw.lock().unwrap().clone(),
                unm.lock().unwrap().clone(),
                diag.lock().unwrap().clone(),
            );
            result
        };
        let single_off = run_single(false);
        let single_on = run_single(true);
        assert_eq!(single_off.0, single_on.0);
        assert_eq!(single_off.1, single_on.1);
        assert_eq!(single_off.2, single_on.2);
        let single_rows: Vec<serde_json::Value> = std::str::from_utf8(&single_on.3)
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        assert_eq!(single_rows.len(), 4);
        for (i, row) in single_rows.iter().enumerate() {
            assert_eq!(row["read_name"], format!("pair{i}/1"));
            assert!(row["mate"].is_null());
            assert!(row["mate_attempt"].is_null());
            assert_eq!(row["pass2_primary_placed"], i == 2);
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    /// pair_records field equality after the borrow-only refactor: flags,
    /// mate coordinates, tlen sign and MAPQ must match the relation
    /// functions exactly (historical construction), for a mapped pair and a
    /// mate-unmapped pair.
    #[test]
    fn pair_records_flags_mates_tlen_match_relation() {
        use crate::extend::CigarOp;
        let aln = |strand: Strand, pos: u32| ReadAlignment {
            contig: 0,
            pos,
            strand,
            score: 90,
            second_score: 0,
            chain_score: 90,
            second_chain_score: 0,
            cigar: vec![CigarOp::Match(50)],
            n_anchors: 10,
            junctions: Vec::new(),
            ea_ag: 0,
            ea_tc: 0,
            mm_count: 0,
            n_seeds: 10,
            rescued: false,
            mapq_override: None,
            seed_search_limited: false,
            exact_path_ambiguous: false,
            short_body_qualified: false,
            pair_substitution_qualified: false,
            origin: crate::mapq::Origin::Genomic,
            n_cands: 1,
        };
        let seq = b"ACGT".repeat(13); // 52
        let qual = vec![b'I'; 52];
        let a1 = aln(Strand::Plus, 100);
        let a2 = aln(Strand::Minus, 400);
        let (ra, rb) = pair_records(
            b"t/1",
            &seq,
            &qual,
            &Some(a1.clone()),
            b"t/2",
            &seq,
            &qual,
            &Some(a2.clone()),
        );
        let proper = relate(&a1, &a2) == PairRelation::Proper;
        let expected_tlen = if proper { template_length(&a1, &a2) } else { 0 };
        assert_eq!(ra.name, b"t/1");
        assert_eq!(rb.name, b"t/2");
        assert!(ra.flag & flag::PAIRED != 0);
        assert!(ra.flag & flag::READ1 != 0);
        assert!(rb.flag & flag::READ2 != 0);
        assert_eq!(ra.flag & flag::PROPER_PAIR != 0, proper);
        assert_eq!(rb.flag & flag::PROPER_PAIR != 0, proper);
        assert_eq!(ra.mate, Some((0, 400, expected_tlen)));
        assert_eq!(rb.mate, Some((0, 100, -expected_tlen)));
        // a2 minus: MATE_REVERSE on R1, REVERSE on R2; R1 plus keeps orientation
        assert!(ra.flag & flag::MATE_REVERSE != 0);
        assert!(rb.flag & flag::REVERSE != 0);
        assert_eq!(ra.flag & flag::REVERSE, 0);
        assert_eq!(rb.seq, esperanto_bamio::revcomp(&seq));
        assert_eq!(ra.seq, seq);
        assert_eq!(ra.mapq, mapq_of(&a1));
        assert_eq!(rb.mapq, mapq_of(&a2));

        // mate unmapped: MATE_UNMAPPED on the partner, no mate coords, tlen 0
        let (rc, rd) = pair_records(b"t/1", &seq, &qual, &Some(a1), b"t/2", &seq, &qual, &None);
        assert!(rc.flag & flag::MATE_UNMAPPED != 0);
        assert_eq!(rc.mate, Some((-1, -1, 0)));
        assert!(rd.flag & flag::UNMAPPED != 0);
        assert_eq!(rd.mapq, 0);
        assert!(rd.aln.is_none());
    }

    /// L1-on must not skip pass-1 junction discovery: the annotation covers
    /// only ~75% of the junctions a sample expresses, so a novel junction is
    /// found by pass 1 and merged into the library for pass 2. Builds a tiny
    /// world (one 20 kb contig, a one-transcript GTF for the L1 bundle, three
    /// pairs whose R1 spans an UNANNOTATED GT..AG intron) and asserts the
    /// discovery is returned in `kept`.
    #[test]
    fn l1_on_still_discovers_novel_junction_in_pass1() {
        // deterministic pseudo-random reference, no repeats
        let n = 20_000usize;
        let mut seq = Vec::with_capacity(n);
        let mut state: u64 = 0x9E3779B97F4A7C15;
        for _ in 0..n {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // novel intron [8000, 9000): GT donor, AG acceptor (plus strand)
        seq[8000] = b'G';
        seq[8001] = b'T';
        seq[8998] = b'A';
        seq[8999] = b'G';
        let dir = std::env::temp_dir().join(format!("esp-2pass-l1-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("ref.fa");
        let mut fasta_txt = b">chrT\n".to_vec();
        for chunk in seq.chunks(60) {
            fasta_txt.extend_from_slice(chunk);
            fasta_txt.push(b'\n');
        }
        std::fs::write(&fa, &fasta_txt).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), b"chrT\t20000\t6\t60\t61\n").unwrap();
        let gtf = dir.join("anno.gtf");
        std::fs::write(
            &gtf,
            b"chrT\tt\ttranscript\t12001\t16000\t.\t+\t.\ttranscript_id \"txA\";\n\
              chrT\tt\texon\t12001\t13000\t.\t+\t.\ttranscript_id \"txA\";\n\
              chrT\tt\texon\t15001\t16000\t.\t+\t.\ttranscript_id \"txA\";\n",
        )
        .unwrap();
        let tidx = dir.join("ref.tidx");
        esperanto_tidx::build(&gtf, &fa, &tidx, &esperanto_tidx::BuildOptions::default()).unwrap();
        let l1 = esperanto_engine::L1Index::build(&tidx, &gtf, &fa).unwrap();

        let reference = crate::fasta::parse_fasta_bytes(&fasta_txt).unwrap();
        let index = crate::index::Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });

        // three PE pairs; R1 spans the novel junction 7925..8000 + 9000..9075
        let mut r1_txt = Vec::new();
        let mut r2_txt = Vec::new();
        for i in 0..3 {
            let mut span = seq[7925..8000].to_vec();
            span.extend_from_slice(&seq[9000..9075]);
            let mate = &seq[9100 + i * 200..9250 + i * 200];
            let q = b"I".repeat(150);
            r1_txt.extend_from_slice(
                &[
                    &format!("@novel{i}/1\n").into_bytes()[..],
                    &span,
                    b"\n+\n",
                    &q,
                    b"\n",
                ]
                .concat(),
            );
            r2_txt.extend_from_slice(
                &[
                    &format!("@novel{i}/2\n").into_bytes()[..],
                    mate,
                    b"\n+\n",
                    &q,
                    b"\n",
                ]
                .concat(),
            );
        }
        let r1 = dir.join("r1.fq");
        let r2 = dir.join("r2.fq");
        std::fs::write(&r1, &r1_txt).unwrap();
        std::fs::write(&r2, &r2_txt).unwrap();

        let mut out = PipelineOut {
            bam: None,
            unmapped_fq: Box::new(std::io::sink()),
            diagnostics: None,
            index: &index,
            config: crate::align::AlignConfig::rna_default(),
            jlib: None,
            l1: Some(std::sync::Arc::new(l1)),
        };
        let (_stats, kept, _merged, _removed) = run_pe_2pass(&mut out, &r1, &r2, 2).unwrap();
        assert!(
            kept.iter().any(|(j, support, fate)| {
                *support >= 2
                    && *fate != DiscoveryFate::OutvotedByAnnotated
                    && !j.minus_strand
                    && (j.start as i64 - 8000).abs() <= 2
                    && (j.end as i64 - 9000).abs() <= 2
            }),
            "novel junction [8000, 9000) must be discovered under L1; kept = {kept:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The pass-1 support vote at an annotated boundary: a shifted
    /// re-observation that out-votes the annotated copy (>2x pass-1 support)
    /// replaces it in the merged library; one that loses is reported but not
    /// merged; within a 2:1 band both stay. Novel junctions (>3 bp from any
    /// annotated one) always stay.
    #[test]
    fn filter_discoveries_support_vote() {
        let annotated = JunctionLib::build([Junction {
            contig: 0,
            start: 100,
            end: 200,
            minus_strand: false,
        }]);
        let j = |start, end| Junction {
            contig: 0,
            start,
            end,
            minus_strand: false,
        };
        // annotated (100,200) observed 40 times in pass 1 (exact re-observation)
        let mut discovered: Discoveries = vec![(j(100, 200), 40, 40 * 60)];
        // shifted +1 with support 7: loses the vote (40 > 2*7)
        discovered.push((j(101, 200), 7, 7 * 60));
        // novel junction: kept
        discovered.push((j(100, 900), 3, 3 * 60));
        // support 1: quality-gated out entirely
        discovered.push((j(5000, 6000), 1, 60));
        let (got, removed) = filter_discoveries(discovered, Some(&annotated));
        assert_eq!(got.len(), 3, "support-1 junction must be gated out");
        assert!(got
            .iter()
            .any(|(j, _, f)| j.start == 100 && j.end == 200 && *f == DiscoveryFate::TieKept));
        assert!(got
            .iter()
            .any(|(j, _, f)| j.start == 101 && *f == DiscoveryFate::OutvotedByAnnotated));
        assert!(got
            .iter()
            .any(|(j, _, f)| j.end == 900 && *f == DiscoveryFate::Kept));
        assert!(removed.is_empty());

        // the shifted copy out-votes an unobserved annotated boundary
        let discovered2: Discoveries = vec![(j(101, 200), 40, 40 * 60)];
        let (got2, removed2) = filter_discoveries(discovered2, Some(&annotated));
        assert!(got2
            .iter()
            .any(|(j, _, f)| j.start == 101 && *f == DiscoveryFate::OutvotesAnnotated));
        assert_eq!(removed2.len(), 1);
        assert_eq!(removed2[0], j(100, 200));

        // without an annotation every survivor is novel-kept
        let (got3, removed3) = filter_discoveries(vec![(j(101, 200), 7, 7 * 60)], None);
        assert!(got3.iter().all(|(_, _, f)| *f == DiscoveryFate::Kept));
        assert!(removed3.is_empty());
    }

    /// Shared unmapped.fq sink: keeps the bytes run_pe wrote so the test can
    /// inspect them after the pipeline finished.
    #[derive(Clone)]
    struct SharedSink(Arc<std::sync::Mutex<Vec<u8>>>);

    impl Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().expect("sink lock").extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// dedup_pe contract: unique table in first-seen order, per-pair uid maps
    /// in chunk order, one uid per sequence across pairs AND across mates.
    #[test]
    fn dedup_pe_uid_mapping_stable_and_order_preserved() {
        let rec = |name: &str, seq: &str| FastqRecord {
            name: name.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: b"I".repeat(seq.len()),
        };
        // pair0/pair2 share the R1 seq; pair0's R1 seq also appears as
        // pair1's R2 seq — the unified table must give it ONE uid.
        let chunk = vec![
            (rec("a/1", "AAAA"), rec("a/2", "CCCC")),
            (rec("b/1", "GGGG"), rec("b/2", "AAAA")),
            (rec("c/1", "AAAA"), rec("c/2", "TTTT")),
        ];
        let (unique_seqs, uid1, uid2) = dedup_pe(&chunk);
        // unique table: first-seen order, each seq exactly once
        assert_eq!(
            unique_seqs,
            vec![&b"AAAA"[..], &b"CCCC"[..], &b"GGGG"[..], &b"TTTT"[..]]
        );
        // one uid per pair position, in chunk order
        assert_eq!(uid1, vec![0, 2, 0]);
        assert_eq!(uid2, vec![1, 0, 3]);
        // same seq → same uid, across pairs (R1 dup) and across mates
        assert_eq!(uid1[0], uid1[2]);
        assert_eq!(uid1[0], uid2[1]);
        // empty chunk: empty tables
        let (u, m1, m2) = dedup_pe(&[]);
        assert!(u.is_empty());
        assert!(m1.is_empty());
        assert!(m2.is_empty());
    }

    /// Paired-end deduplication must not change per-pair outputs: a chunk whose
    /// sequences are duplicated (a second block of pairs with identical
    /// R1/R2 seqs but different names and quals) must yield exactly the same
    /// per-pair results as the unduplicated chunk — every stats counter
    /// doubles, the discovery stream is the exact block concatenation, and
    /// the unmapped.fq seq multiset doubles (quality must not move the
    /// primary alignment).
    #[test]
    fn pe_dedup_duplicated_chunk_matches_unduplicated_outputs() {
        // deterministic pseudo-random reference, no repeats (same recipe as
        // the L1 discovery test)
        let n = 20_000usize;
        let mut seq = Vec::with_capacity(n);
        let mut state: u64 = 0x9E3779B97F4A7C15;
        for _ in 0..n {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // GT..AG intron [8000, 9000) so one pair exercises the junction path
        seq[8000] = b'G';
        seq[8001] = b'T';
        seq[8998] = b'A';
        seq[8999] = b'G';
        // an R2 with no k-mer in the reference (LCG continuation — fixed
        // bytes, so non-matching is a fixed fact, not a flake)
        let mut absent = Vec::with_capacity(150);
        for _ in 0..150 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            absent.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // pairs: (R1, R2). p1's R1 duplicates p0's R1; p2's short R1 makes
        // its own placement poor (rescue attempt consuming the cached
        // placements); p4's absent R2 gives the unplaced-mate rescue path
        // (cache-free).
        let spliced = {
            let mut s = seq[7925..8000].to_vec();
            s.extend_from_slice(&seq[9000..9075]);
            s
        };
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = vec![
            (seq[1000..1150].to_vec(), seq[1200..1350].to_vec()),
            (seq[1000..1150].to_vec(), seq[3000..3150].to_vec()),
            (seq[5000..5080].to_vec(), seq[4000..4150].to_vec()),
            (spliced, seq[9100..9250].to_vec()),
            (seq[7000..7150].to_vec(), absent),
        ];
        let dir = std::env::temp_dir().join(format!("esp-pe-dedup-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let rec_bytes = |name: &str, s: &[u8], qual: u8| -> Vec<u8> {
            [
                &format!("@{name}\n").into_bytes()[..],
                s,
                b"\n+\n",
                &vec![qual; s.len()],
                b"\n",
            ]
            .concat()
        };
        let mut base1 = Vec::new();
        let mut base2 = Vec::new();
        let mut dup1 = Vec::new();
        let mut dup2 = Vec::new();
        for (i, (s1, s2)) in pairs.iter().enumerate() {
            base1.extend(rec_bytes(&format!("p{i}"), s1, b'I'));
            base2.extend(rec_bytes(&format!("p{i}"), s2, b'I'));
            dup1.extend(rec_bytes(&format!("p{i}"), s1, b'I'));
            dup2.extend(rec_bytes(&format!("p{i}"), s2, b'I'));
        }
        // duplicated block: same seqs, different names and quals
        for (i, (s1, s2)) in pairs.iter().enumerate() {
            dup1.extend(rec_bytes(&format!("q{i}"), s1, b'J'));
            dup2.extend(rec_bytes(&format!("q{i}"), s2, b'J'));
        }
        let base_r1 = dir.join("base_1.fq");
        let base_r2 = dir.join("base_2.fq");
        let dup_r1 = dir.join("dup_1.fq");
        let dup_r2 = dir.join("dup_2.fq");
        std::fs::write(&base_r1, &base1).unwrap();
        std::fs::write(&base_r2, &base2).unwrap();
        std::fs::write(&dup_r1, &dup1).unwrap();
        std::fs::write(&dup_r2, &dup2).unwrap();

        let reference = crate::fasta::parse_fasta_bytes(&{
            let mut t = b">chrT\n".to_vec();
            for chunk in seq.chunks(60) {
                t.extend_from_slice(chunk);
                t.push(b'\n');
            }
            t
        })
        .unwrap();
        let index = crate::index::Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });

        let run = |r1: &Path, r2: &Path| -> (AlignStats, Discoveries, Vec<u8>) {
            let cell = Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut out = PipelineOut {
                bam: None,
                unmapped_fq: Box::new(SharedSink(cell.clone())),
                diagnostics: None,
                index: &index,
                config: crate::align::AlignConfig::rna_default(),
                jlib: None,
                l1: None,
            };
            let (stats, disc) = run_pe(&mut out, r1, r2, 2).unwrap();
            let raw = std::mem::take(&mut *cell.lock().unwrap());
            (stats, disc, raw)
        };
        let (s1, d1, raw1) = run(&base_r1, &base_r2);
        let (s2, d2, raw2) = run(&dup_r1, &dup_r2);

        // the base run is non-trivial: something maps, the spliced pair
        // discovers the junction, the absent-R2 pair leaves unmapped output
        assert!(s1.mapped_reads > 0, "base run must map something");
        assert!(!d1.is_empty(), "spliced pair must discover the junction");
        let unm_seq_lines = |raw: &[u8]| -> Vec<Vec<u8>> {
            use std::io::Read as _;
            let mut txt = Vec::new();
            flate2::read::GzDecoder::new(raw)
                .read_to_end(&mut txt)
                .unwrap();
            txt.split(|&b| b == b'\n')
                .enumerate()
                .filter(|(i, _)| i % 4 == 1)
                .map(|(_, l)| l.to_vec())
                .collect()
        };
        let mut u1 = unm_seq_lines(&raw1);
        assert!(!u1.is_empty(), "absent R2 must stay unmapped");

        // every observable output doubles exactly
        assert_eq!(s2.total_reads, 2 * s1.total_reads);
        assert_eq!(s2.mapped_reads, 2 * s1.mapped_reads);
        assert_eq!(s2.unmapped_reads, 2 * s1.unmapped_reads);
        assert_eq!(s2.proper_pairs, s1.proper_pairs.map(|x| 2 * x));
        assert_eq!(s2.junctions_total, 2 * s1.junctions_total);
        assert_eq!(s2.rescued_total, 2 * s1.rescued_total);
        assert_eq!(s2.rescue_fail_total, 2 * s1.rescue_fail_total);
        assert_eq!(s2.decoy_reads, 2 * s1.decoy_reads);
        for (a, b) in s1.mapq_hist.iter().zip(s2.mapq_hist.iter()) {
            assert_eq!(2 * a, *b, "mapq histogram must double");
        }
        // discovery stream: exact block duplication (pair order preserved)
        assert_eq!(d2.len(), 2 * d1.len());
        assert_eq!(&d2[..d1.len()], &d1[..]);
        assert_eq!(&d2[d1.len()..], &d1[..]);
        // unmapped.fq: same seq multiset, doubled
        let mut u2 = unm_seq_lines(&raw2);
        u1.sort();
        u2.sort();
        let mut doubled = u1.clone();
        doubled.extend_from_slice(&u1);
        assert_eq!(u2, doubled);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
