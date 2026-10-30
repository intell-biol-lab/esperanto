//! Per-read alignment driver: seed → anchor → chain → (L1 / intron-chain /
//! spliced / EA-Myers / banded DP) → split/span rescue → recount → coverage gate.
//!
//! Alignment gates include: the splice gate (skip try_spliced on
//! <2 segments), the non-best-chain spliced Myers pre-gate (skip the spliced DP when the
//! pseudo-reference is farther than read_len/3 edits from the read), the tail-rescue gate
//! (small-tail rescue tightening), multitie marking, EA-Myers decisive path, window cap.

use std::collections::HashSet;
use std::sync::Arc;

use crate::chain::{
    chain_anchors_into, second_score_with, tied_best, Chain, ChainScratch, TOP_CANDIDATES,
};
use crate::extend::{
    extend_hint, CigarOp, DiagHint, EditDir, ExtendBuffer, ExtendParams, Extension,
};
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::index::{collect_anchors_into, Index};
use crate::mapq::ReadAlignment;
use crate::seed::{
    append_query_edge_kmers, minimizers_into, Minimizer, SeedParams, SeedScratch, Strand,
};
use crate::splice::{
    align_spliced, intron_length_penalty, segment_chain, splice_signal, SpliceParams,
};

/// Keep uncovered query ends available without increasing the maximum
/// continuous-reference window width. Both verification paths use this geometry.
fn continuous_extension_window(
    chain: &Chain,
    k: u32,
    read_len: u32,
    contig_len: u32,
) -> (u32, u32) {
    let flank = read_len / 2 + 30;
    let (rs, re) = chain.ref_span(k);
    let (qs, qe) = chain.read_span(k);
    let head = flank.max(qs);
    let tail = flank.max(read_len.saturating_sub(qe));
    let start = rs.saturating_sub(head).min(contig_len);
    let end = re
        .saturating_add(tail)
        .min(contig_len)
        .min(rs.saturating_add(read_len).saturating_add(flank))
        .min(
            start
                .saturating_add(read_len)
                .saturating_add(flank.saturating_mul(2)),
        );
    (start, end.max(start))
}

/// Aligner configuration (`default` = DNA preset; `rna_default` = RNA production).
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    pub seed: SeedParams,
    pub chain: crate::chain::ChainParams,
    pub extend: ExtendParams,
    /// RNA mode: dense seeding + intron-aware chaining + splice machinery.
    pub rna: bool,
    /// Infer conservative paired read-through endpoints before sequence caching.
    pub fragment_endpoints: bool,
    /// Additional bounded high-frequency ordinary context hypotheses.
    pub repeat_context: bool,
    /// Require additional context proposals to beat a constant-base background.
    pub context_background: bool,
    /// Experimental full-query exact-path arbitration (disabled by default).
    pub exact_path_competition: bool,
    /// Qualify otherwise unplaced exact short bodies with observable pair evidence.
    pub mate_short_body: bool,
}

impl Default for AlignConfig {
    /// DNA preset (k=15, w=10).
    fn default() -> Self {
        AlignConfig {
            seed: SeedParams { k: 15, w: 10 },
            chain: crate::chain::ChainParams::default(),
            extend: ExtendParams::default(),
            rna: false,
            fragment_endpoints: false,
            repeat_context: false,
            context_background: false,
            exact_path_competition: false,
            mate_short_body: false,
        }
    }
}

impl AlignConfig {
    /// RNA mode defaults (spec §parameter table).
    pub fn rna_default() -> Self {
        AlignConfig {
            seed: SeedParams { k: 15, w: 5 },
            chain: crate::chain::ChainParams::rna_default(),
            extend: ExtendParams::default(),
            rna: true,
            fragment_endpoints: true,
            repeat_context: true,
            context_background: true,
            exact_path_competition: false,
            mate_short_body: true,
        }
    }
}

/// Second-round (mid-occ) seed occurrence cap (minimap2 --mid-occ semantics).
pub const MID_OCC_CAP: u32 = 5000;
/// Max candidates entering base-level tie-breaking (deterministic truncation).
pub const MAX_TIED_CANDIDATES: usize = 8;
/// Extension window hard cap (garbage chains must not blow up the DP).
pub const MAX_EXTEND_WINDOW: u32 = 20_000;
/// Extension tie-break budget (top-N candidates extended by banded DP).
pub const EXTEND_BUDGET: usize = 2;

/// Outcome of one seeding+chaining round.
pub(crate) struct RoundOutcome {
    chains: Vec<Chain>,
    second: i32,
}

impl RoundOutcome {
    fn best(&self) -> Option<&Chain> {
        self.chains.first()
    }
}

/// Minimum fraction of a read's bases a placement must explain, below which
/// Tier-1 refuses it outright (STAR's `outFilterMatchNminOverLread`
/// community default). Refused reads are unmapped and feed the Tier-2
/// collapsed-alphabet rescue pool.
pub const REFUSE_MIN_ALIGNED_FRAC: f64 = 0.66;

/// Read-span coverage of the round-1 best chain below which the
/// mid-occurrence seeding round runs (see `genomic_candidates`).
pub const MID_ROUND_MIN_COVERAGE: f64 = 0.6;

/// Anchor budget for the mid-occurrence round. On overflow, admit whole
/// seed hit lists, rarest first, and mark the search incomplete (MAPQ 0).
pub const MID_ROUND_MAX_ANCHORS: usize = 20_000;

/// Half-width of the mate-anchored window, on top of the insert estimate
/// (`est_insert + est_insert/2 + read_len` each side of the mate).
pub const MATE_WINDOW_INSERT: u32 = 500;

/// Strict-domain score of one placement candidate (see
/// [`Aligner::rescore_strict`]). All candidates of a read are scored by the
/// same function with the same parameters, so these numbers are comparable
/// across layers and decide both placement and MAPQ.
#[derive(Clone, Debug, Default)]
pub(crate) struct StrictScore {
    /// Total strict score (match/mismatch/indel/splice-signal sums).
    pub score: i32,
    /// Read bases in M/I columns.
    pub aligned: u32,
    /// Exact matches in M columns.
    pub matches: u32,
    /// Ordinary mismatches (N on either side excluded).
    pub mm: u32,
    /// A-to-G-direction substitutions (ref A / read G), charged as mismatches.
    pub ea_ag: u32,
    /// T-to-C-direction substitutions (ref T / read C), charged as mismatches.
    pub ea_tc: u32,
    /// Insertion events.
    pub ins: u32,
    /// Deletion events.
    pub del: u32,
    /// Soft-clipped bases.
    pub clipped: u32,
    /// Intron (RefSkip) count.
    pub n_splice: u32,
    /// Junction prior contribution, separate from base/indel evidence.
    pub splice_score: i32,
}

/// One placement candidate from either layer, pending strict arbitration.
#[derive(Clone)]
pub(crate) struct Placement {
    pub origin: crate::mapq::Origin,
    pub contig: u32,
    pub pos: u32,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub junctions: Vec<RefinedJunction>,
    /// Chain score of the source chain (genomic layer; 0 for L1).
    pub chain_score: i32,
    /// Score the producing extension reported (domain: extension scores).
    pub ext_score: i32,
    pub strict: StrictScore,
    /// Editing direction the candidate was last scored under (`Off` = strict).
    pub ea_dir: EditDir,
}

impl Placement {
    /// Query evidence, not reference span: strict.aligned historically also
    /// includes deleted reference columns. Preserve that score bookkeeping
    /// for MAPQ/edit-density consumers, but never credit D/N as read coverage.
    pub(crate) fn informative_query_coverage(&self) -> u32 {
        let deleted: u32 = self
            .cigar
            .iter()
            .filter_map(|op| match op {
                CigarOp::Del(n) => Some(*n),
                _ => None,
            })
            .sum();
        self.strict.aligned.saturating_sub(deleted)
    }

    fn coverage_eligible(&self, read_len: usize) -> bool {
        let covered = self.informative_query_coverage();
        read_len > 0 && f64::from(covered) >= REFUSE_MIN_ALIGNED_FRAC * read_len as f64
    }

    pub(crate) fn primary_eligible(&self, read_len: usize) -> bool {
        self.strict.score > 0 && self.coverage_eligible(read_len)
    }

    /// A fast proposal's skipped intervals must be the junctions it reports.
    /// Changing N length alone cannot move its donor or read breakpoint.
    fn has_consistent_junctions(&self) -> bool {
        let mut observed = Vec::new();
        let mut rf = self.pos;
        for op in &self.cigar {
            match op {
                CigarOp::Match(n) | CigarOp::Del(n) => rf += *n,
                CigarOp::RefSkip(n) => {
                    observed.push((self.contig, rf, rf + *n));
                    rf += *n;
                }
                _ => {}
            }
        }
        let mut declared: Vec<_> = self
            .junctions
            .iter()
            .map(|j| (j.junction.contig, j.junction.start, j.junction.end))
            .collect();
        observed.sort_unstable();
        declared.sort_unstable();
        observed == declared
    }

    fn from_genomic(
        contig: u32,
        pos: u32,
        strand: Strand,
        cigar: Vec<CigarOp>,
        junctions: Vec<RefinedJunction>,
        chain_score: i32,
        ext_score: i32,
    ) -> Placement {
        Placement {
            origin: crate::mapq::Origin::Genomic,
            contig,
            pos,
            strand,
            cigar,
            junctions,
            chain_score,
            ext_score,
            strict: StrictScore::default(),
            ea_dir: EditDir::Off,
        }
    }
}

/** Genomic layer output: per-chain placement candidates plus read evidence.
 */
pub(crate) struct GenomicBatch {
    pub placements: Vec<Placement>,
    pub n_seeds: usize,
    pub seed_budget: Option<crate::index::AnchorBudget>,
}

/** Per-thread aligner state (reuses DP buffers).
 */
pub struct Aligner<'a> {
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb / 2-pass discoveries; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// L1 transcriptome-first engine bundle (None = genomic alignment only).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
    pub(crate) buf: ExtendBuffer,
    /// Reusable hot-path buffers (seeding / chaining / revcomp / windows).
    /// Allocation-elimination only: every buffer is cleared before refill,
    /// so results are byte-identical to the fresh-allocation path.
    pub(crate) scratch: MapScratch,
}

/** Reusable buffers for the seeding/chaining hot path. Owned by the
 * [`Aligner`] (one per rayon worker, via `map_init`) so per-call `Vec`
 * allocations become `clear()` + refill. No algorithm state lives here.
 */
#[derive(Default)]
pub(crate) struct MapScratch {
    seed: SeedScratch,
    chain: ChainScratch,
    second_keys: HashSet<(u32, u32)>,
    anchors: Vec<crate::seed::Anchor>,
    mins_out: Vec<Minimizer>,
    chains_out: Vec<Chain>,
    rc: Vec<u8>,
    win: Vec<u8>,
    myers_a: Vec<u8>,
    myers_b: Vec<u8>,
    rescore: Vec<u8>,
    kmer_scan: Vec<u8>,
}

/** The scored candidate list of one `align_read` call, kept for the
 * mate-anchored re-arbitration so the read is not aligned twice.
 */
#[derive(Clone)]
pub struct CachedPlacements {
    placements: Vec<Placement>,
    n_seeds: usize,
    genomic_any: bool,
    seed_budget: Option<crate::index::AnchorBudget>,
}

impl CachedPlacements {
    pub(crate) fn seed_search_limited(&self) -> bool {
        self.seed_budget.is_some_and(|b| b.limited())
    }

    pub(crate) fn diagnostic_pool(
        &self,
        index: &Index,
        read_len: usize,
        stage: &'static str,
    ) -> crate::diagnostics::CandidatePool {
        let mut pool =
            crate::diagnostics::CandidatePool::capture(index, read_len, &self.placements, stage);
        pool.seed_budget = self.seed_budget;
        pool
    }
}

impl<'a> Aligner<'a> {
    pub fn new(index: &'a Index, config: AlignConfig) -> Self {
        Aligner {
            index,
            config,
            jlib: None,
            l1: None,
            buf: ExtendBuffer::default(),
            scratch: MapScratch::default(),
        }
    }

    pub fn with_lib(mut self, lib: Arc<JunctionLib>) -> Self {
        self.jlib = Some(lib);
        self
    }

    /// Enable the L1 transcriptome-first fast path: hit → direct placement;
    /// miss (Fallback) → genomic alignment.
    pub fn with_l1(mut self, l1: Arc<esperanto_engine::L1Index>) -> Self {
        self.l1 = Some(l1);
        self
    }

    /// Take the reusable reverse-complement buffer, filled with
    /// `revcomp(seq)` (allocation-elimination only; contents identical).
    /// Give it back with `self.scratch.rc = v` when done.
    fn take_rc(&mut self, seq: &[u8]) -> Vec<u8> {
        let mut v = std::mem::take(&mut self.scratch.rc);
        revcomp_into(seq, &mut v);
        v
    }

    /// Take the reusable reference-window buffer, decoded over
    /// `[wstart, wend)` (identical bytes to `contig.slice_ascii`).
    fn take_win(&mut self, wstart: u32, wend: u32, contig: &crate::fasta::Contig) -> Vec<u8> {
        let mut v = std::mem::take(&mut self.scratch.win);
        v.clear();
        contig.decode_into(wstart, wend, &mut v);
        v
    }

    /// Convert one L1 engine candidate into a genomic-layer placement.
    ///
    /// Conservation guard: the projected CIGAR's query consumption (M+I+S)
    /// must equal read_len, otherwise the projection is discarded.
    fn l1_placement(&self, read_len: usize, c: esperanto_engine::L1Candidate) -> Option<Placement> {
        let genome_contig = self
            .l1
            .as_ref()
            .and_then(|l1| l1.contig_name(c.contig))
            .and_then(|name| self.index.reference.contig_index(name.as_bytes()))?;

        let strand = match c.strand {
            esperanto_engine::Strand::Plus => Strand::Plus,
            esperanto_engine::Strand::Minus => Strand::Minus,
        };
        let cigar: Vec<CigarOp> = c
            .cigar
            .iter()
            .map(|op| match op {
                esperanto_engine::CigarOp::Match(n) => CigarOp::Match(*n),
                esperanto_engine::CigarOp::Ins(n) => CigarOp::Ins(*n),
                esperanto_engine::CigarOp::Del(n) => CigarOp::Del(*n),
                esperanto_engine::CigarOp::RefSkip(n) => CigarOp::RefSkip(*n),
                esperanto_engine::CigarOp::SoftClip(n) => CigarOp::SoftClip(*n),
            })
            .collect();
        let qspan: usize = cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n as usize,
                _ => 0,
            })
            .sum();
        if qspan != read_len {
            return None; // non-conserving projection
        }
        Some(Placement {
            origin: crate::mapq::Origin::L1,
            contig: genome_contig,
            pos: c.pos,
            strand,
            cigar,
            junctions: Vec::new(),
            chain_score: 0,
            ext_score: 0,
            strict: StrictScore::default(),
            ea_dir: EditDir::Off,
        })
    }

    /// Align one read (ASCII ACGTN).
    ///
    /// Tier-1 unified arbitration: both layers act as candidate generators —
    /// the L1 engine contributes every defensible transcriptome projection,
    /// the genomic layer one extension per chained locus — and every
    /// candidate is rescored by one strict four-letter function with the
    /// editing exemption off. The best eligible strict score decides the placement
    /// (ties prefer the genomic layer, whose CIGAR was computed against the
    /// genome directly), the margin to the best non-overlapping runner-up
    /// decides MAPQ. Primary eligibility requires positive score and at least
    /// [`REFUSE_MIN_ALIGNED_FRAC`] informative query coverage; if no candidate
    /// qualifies the read reaches the Tier-2 pool. Ineligible different-locus
    /// evidence still contributes to the ambiguity margin.
    pub fn align_read(&mut self, seq: &[u8]) -> Option<ReadAlignment> {
        let (aln, _) = self.align_read_full(seq);
        aln
    }

    /// `align_read` plus the scored candidate list, for the mate-anchored
    /// re-arbitration to reuse instead of re-aligning the read from scratch.
    pub fn align_read_full(
        &mut self,
        seq: &[u8],
    ) -> (Option<ReadAlignment>, Option<CachedPlacements>) {
        let Some((placements, n_seeds, genomic_any, seed_budget)) = self.scored_placements(seq)
        else {
            return (None, None);
        };
        let cache = CachedPlacements {
            placements: placements.clone(),
            n_seeds,
            genomic_any,
            seed_budget,
        };
        (
            self.arbitrate_with_budget(seq, placements, n_seeds, genomic_any, seed_budget),
            Some(cache),
        )
    }

    /// Collect, strictly rescore, coverage-class dedup and order the placements of
    /// a read from both layers. Retains an empty pool when an incomplete seed
    /// search must be propagated to later mate-window arbitration.
    fn scored_placements(
        &mut self,
        seq: &[u8],
    ) -> Option<(
        Vec<Placement>,
        usize,
        bool,
        Option<crate::index::AnchorBudget>,
    )> {
        let mut placements: Vec<Placement> = Vec::new();

        // 1. Candidate generation — transcriptome layer.
        if let Some(l1) = self.l1.clone() {
            let mut stats = esperanto_engine::ReadStats::default();
            for c in esperanto_engine::align_read_candidates(
                seq,
                l1.as_ref(),
                l1.as_ref(),
                l1.as_ref(),
                &esperanto_engine::EngineConfig::default(),
                &mut stats,
            ) {
                if let Some(p) = self.l1_placement(seq.len(), c) {
                    placements.push(p);
                }
            }
        }

        // 2. Candidate generation — genomic layer (always consulted).
        let batch = self.genomic_candidates(seq);
        let n_seeds = batch.n_seeds;
        let genomic_any = !batch.placements.is_empty();
        let seed_budget = batch.seed_budget;
        placements.extend(batch.placements);
        if placements.is_empty() && seed_budget.is_none() {
            return None;
        }

        // 3. Unified strict rescoring: same function, same parameters, no
        // editing exemption.
        let rc_storage = self.take_rc(seq);
        for p in &mut placements {
            let q: &[u8] = if p.strand == Strand::Minus {
                &rc_storage
            } else {
                seq
            };
            self.rescore_strict(q, p);
        }

        // 4. Same-locus dedup within coverage classes, then decision order.
        // Retain full nonpositive candidates for the hyperediting check below.
        dedup_by_footprint_with_options(
            &mut placements,
            seq.len(),
            self.config.exact_path_competition,
        );
        sort_placements(&mut placements);

        // 5. Hyperedited read: re-arbitrate in its own reduced alphabet (see
        // `rescore_with`).
        self.hyper_rearbitrate(seq, &rc_storage, &mut placements);
        self.scratch.rc = rc_storage;
        Some((placements, n_seeds, genomic_any, seed_budget))
    }

    /// Hyperedited-read re-arbitration (see `rescore_with`). The signature is
    /// looked for in every candidate, not only the strict winner: the strict
    /// extension soft-clips an edited tail (each edit costs 4), so the winner
    /// is often the clipped version whose remaining edit density no longer
    /// qualifies, while a full-length candidate of the same locus (L1
    /// projection, or the clipped one re-extended) shows the signature
    /// plainly. When found, every unspliced genomic candidate is re-extended
    /// and every candidate rescored in that one read-oriented conversion
    /// domain. A candidate on the opposite alignment strand must complement
    /// the reference-oriented direction (A>G <-> T>C). Then dedup and order
    /// are recomputed. No-op otherwise.
    fn hyper_rearbitrate(
        &mut self,
        seq: &[u8],
        rc_storage: &[u8],
        placements: &mut Vec<Placement>,
    ) {
        let Some((dir, evidence_strand)) = placements
            .iter()
            .find_map(|p| hyper_direction(p).map(|dir| (dir, p.strand)))
        else {
            return;
        };
        for p in placements.iter_mut() {
            let candidate_dir = if p.strand == evidence_strand {
                dir
            } else {
                match dir {
                    EditDir::Ag => EditDir::Tc,
                    EditDir::Tc => EditDir::Ag,
                    EditDir::Off => EditDir::Off,
                }
            };
            let q: &[u8] = if p.strand == Strand::Minus {
                rc_storage
            } else {
                seq
            };
            if p.origin == crate::mapq::Origin::Genomic {
                self.reextend_with(q, p, candidate_dir);
            }
            self.rescore_with(q, p, candidate_dir);
        }
        dedup_by_footprint_with_options(placements, seq.len(), self.config.exact_path_competition);
        sort_placements(placements);
    }

    /// Align a read whose independent placement is poor, adding the locus
    /// the mate points at as one more candidate. The mate's window is
    /// searched directly (bit-parallel infix search, then strict extension),
    /// so a read whose every k-mer sits above the occurrence caps still
    /// reaches the locus its fragment belongs to — the way STAR and BWA
    /// resolve a repeat copy by the pair. The candidate is arbitrated with
    /// the read's own candidates under the same strict score; when it wins,
    /// the placement is marked rescued (MAPQ capped at the pileup gate, `RE`
    /// tag), because it was won on the mate's evidence rather than the read's
    /// own seeds. Returns None when the mate-window candidate does not win.
    pub fn align_read_with_mate(
        &mut self,
        seq: &[u8],
        mate_contig: u32,
        mate_pos: u32,
        cached: Option<&CachedPlacements>,
    ) -> Option<ReadAlignment> {
        self.align_read_with_mate_observed(seq, mate_contig, mate_pos, cached, None)
    }

    pub(crate) fn align_read_with_mate_observed(
        &mut self,
        seq: &[u8],
        mate_contig: u32,
        mate_pos: u32,
        cached: Option<&CachedPlacements>,
        mut diagnostic: Option<&mut Option<crate::diagnostics::CandidatePool>>,
    ) -> Option<ReadAlignment> {
        if let Some(slot) = diagnostic.as_mut() {
            **slot = Some(crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &[],
                "mate_window_no_proposal",
            ));
            if let Some(pool) = slot.as_mut() {
                pool.seed_budget = cached.and_then(|c| c.seed_budget);
            }
        }
        let mut win = self.mate_window_placement(seq, mate_contig, mate_pos)?;
        let rc_storage = self.take_rc(seq);
        let q: &[u8] = if win.strand == Strand::Minus {
            &rc_storage
        } else {
            seq
        };
        self.rescore_strict(q, &mut win);
        let win_probe = win.clone();
        let (mut placements, n_seeds, _, seed_budget) = match cached {
            // Reuse the caller's own scored candidates: the read was just
            // aligned once already.
            Some(c) => (
                c.placements.clone(),
                c.n_seeds,
                c.genomic_any,
                c.seed_budget,
            ),
            None => self
                .scored_placements(seq)
                .unwrap_or((Vec::new(), 0, false, None)),
        };
        placements.push(win);
        dedup_by_footprint_with_options(
            &mut placements,
            seq.len(),
            self.config.exact_path_competition,
        );
        self.hyper_rearbitrate(seq, &rc_storage, &mut placements);
        if let Some(slot) = diagnostic.as_mut() {
            **slot = Some(crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "mate_window_rejected",
            ));
            if let Some(pool) = slot.as_mut() {
                pool.seed_budget = seed_budget;
            }
        }
        // The mate-window candidate must have won an eligible locus (dedup keeps the
        // better-scored of two same-locus candidates, so compare footprints).
        if !placements
            .iter()
            .find(|p| p.primary_eligible(seq.len()))
            .is_some_and(|p| footprint_overlap(p, &win_probe) >= 0.5)
        {
            self.scratch.rc = rc_storage;
            return None;
        }
        self.scratch.rc = rc_storage;
        // Pair consistency prunes the competition: a runner-up on another
        // contig, or outside the fragment's reach of the mate, is not a
        // placement of this fragment. MAPQ is then the margin over the
        // mate-consistent alternatives only, which is how a paired aligner
        // resolves a repeat copy by its pair (STAR, BWA-MEM).
        let reach = 2 * (MATE_WINDOW_INSERT + MATE_WINDOW_INSERT / 2 + seq.len() as u32);
        placements.retain(|p| {
            (p.contig == mate_contig && p.pos.abs_diff(mate_pos) <= reach)
                || (self.config.exact_path_competition
                    && footprint_overlap(p, &win_probe) >= 0.5
                    && exact_path::distinct(p, &win_probe, seq.len()))
        });
        if let Some(slot) = diagnostic {
            *slot = Some(crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "mate_window_arbitration",
            ));
            if let Some(pool) = slot.as_mut() {
                pool.seed_budget = seed_budget;
            }
        }
        // The genomic layer was consulted (the window search is genomic), so
        // the L1-only last-resort policy does not apply here.
        let mut aln = self.arbitrate_with_budget(seq, placements, n_seeds, true, seed_budget)?;
        aln.rescued = true;
        Some(aln)
    }

    /// Candidate at the mate's locus: the read (both orientations) is placed
    /// in a window around the mate by the EA-Myers infix search and extended
    /// strictly on the diagonal it finds. No seeds are consulted, so the
    /// occurrence caps cannot hide the locus. None when the window holds no
    /// placement within the distance gate.
    fn mate_window_placement(
        &mut self,
        seq: &[u8],
        mate_contig: u32,
        mate_pos: u32,
    ) -> Option<Placement> {
        let read_len = seq.len() as u32;
        if read_len == 0 || read_len > 256 {
            return None;
        }
        let m = read_len as usize;
        let contig = &self.index.reference.contigs[mate_contig as usize];
        let span = MATE_WINDOW_INSERT + MATE_WINDOW_INSERT / 2 + read_len;
        let wstart = mate_pos.saturating_sub(span);
        let wend = mate_pos.saturating_add(span + read_len).min(contig.len);
        if wend <= wstart || ((wend - wstart) as usize) < m {
            return None;
        }
        let window = self.take_win(wstart, wend, contig);
        let rc_storage = self.take_rc(seq);
        let (mut myers_a, mut myers_b) = (
            std::mem::take(&mut self.scratch.myers_a),
            std::mem::take(&mut self.scratch.myers_b),
        );
        let threshold = ((read_len as i32) / 3).max(10);
        let mut best: Option<(Strand, i32, usize)> = None;
        for (strand, q) in [(Strand::Plus, seq), (Strand::Minus, &rc_storage[..])] {
            let (d, start) = if m <= 128 {
                crate::myers_ea::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
            } else {
                crate::myers_ea::long::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
            };
            if d > threshold || start + m > window.len() {
                continue;
            }
            if best.as_ref().is_none_or(|b| d < b.1) {
                best = Some((strand, d, start));
            }
        }
        self.scratch.myers_a = myers_a;
        self.scratch.myers_b = myers_b;
        let Some((strand, _, start)) = best else {
            self.scratch.rc = rc_storage;
            self.scratch.win = window;
            return None;
        };
        let q: &[u8] = match strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let hint = DiagHint {
            offset: start as i64,
            num: 1,
            den: 1,
        };
        let ext = extend_hint(q, &window, &self.config.extend, &mut self.buf, hint);
        self.scratch.rc = rc_storage;
        self.scratch.win = window;
        if ext.score <= 0 {
            return None;
        }
        Some(Placement::from_genomic(
            mate_contig,
            wstart + ext.ref_start,
            strand,
            ext.cigar,
            Vec::new(),
            0,
            ext.score,
        ))
    }

    /// The Tier-1 decision over an ordered, deduplicated candidate list.
    fn arbitrate_with_budget(
        &mut self,
        seq: &[u8],
        placements: Vec<Placement>,
        n_seeds: usize,
        genomic_any: bool,
        seed_budget: Option<crate::index::AnchorBudget>,
    ) -> Option<ReadAlignment> {
        let mut aln = self.arbitrate(seq, placements, n_seeds, genomic_any)?;
        if seed_budget.is_some_and(|b| b.limited()) {
            aln.seed_search_limited = true;
            crate::mapq::force_mapq0(&mut aln);
        }
        Some(aln)
    }

    fn arbitrate(
        &mut self,
        seq: &[u8],
        placements: Vec<Placement>,
        n_seeds: usize,
        genomic_any: bool,
    ) -> Option<ReadAlignment> {
        // A score leader that lacks read coverage cannot veto a lower-ranked
        // eligible placement. Keep it as ambiguity evidence at other loci.
        let best = placements.iter().find(|p| p.primary_eligible(seq.len()))?;
        // Annotation priors cannot distinguish two complete exact mappings
        // of this query. Keep their actual scores, but display a deterministic
        // representative and preserve the unresolved path evidence.
        let (best, exact_path_ambiguous) = if self.config.exact_path_competition {
            exact_path::representative(&placements, best, seq.len())
        } else {
            (best, false)
        };
        // Runner-up must be a different locus (footprint overlap < 0.5).
        let second = placements
            .iter()
            .filter(|p| footprint_overlap(best, p) < 0.5)
            .map(|p| p.strict.score)
            .filter(|s| *s > 0)
            .max()
            .unwrap_or(0);

        // A best hit on a decoy contig is the decoy contract doing its job:
        // the read's source is absent from the primary reference (rDNA,
        // retrotransposed copies, contamination), and parking it on the decoy
        // keeps it out of the pileup while remaining visible and auditable in
        // the BAM. Decoy placements never enter the Tier-2 pool.
        let best_origin = if is_decoy_contig(self.index.reference, best.contig) {
            crate::mapq::Origin::Decoy
        } else {
            best.origin
        };

        let rc_storage = self.take_rc(seq);
        let q: &[u8] = if best.strand == Strand::Minus {
            &rc_storage
        } else {
            seq
        };
        let n_anchors = self.exact_kmer_support(q, best);

        let mut aln = ReadAlignment {
            contig: best.contig,
            pos: best.pos,
            strand: best.strand,
            score: best.strict.score,
            second_score: second,
            chain_score: if best.origin == crate::mapq::Origin::L1 {
                best.strict.score
            } else {
                best.chain_score
            },
            second_chain_score: 0,
            cigar: best.cigar.clone(),
            n_anchors,
            junctions: best.junctions.clone(),
            ea_ag: best.strict.ea_ag,
            ea_tc: best.strict.ea_tc,
            mm_count: best.strict.mm,
            n_seeds,
            rescued: false,
            mapq_override: None,
            origin: best_origin,
            n_cands: placements
                .iter()
                .filter(|p| p.primary_eligible(seq.len()))
                .count(),
            seed_search_limited: false,
            exact_path_ambiguous,
            short_body_qualified: false,
        };

        // Only the transcriptome layer could place the read. Its apparent
        // ambiguity has a cause when the read carries a hyperediting
        // signature — dense single-direction A>G (or T>C) substitution
        // destroys the exact-seed discrimination the genomic layer runs on,
        // which is the whole reason the Tier-2 rescue exists. Such a
        // placement takes the pileup-gate ceiling; anything else the
        // transcriptome layer alone placed is a last resort and scores 0.
        if !genomic_any && aln.origin == crate::mapq::Origin::L1 {
            if aln.hyperedited_signature() {
                aln.mapq_override = Some(crate::mapq::PILEUP_GATE);
            } else {
                crate::mapq::force_mapq0(&mut aln);
            }
        }
        self.scratch.rc = rc_storage;

        Some(aln)
    }

    /// Rescore one candidate in the strict domain: per-base match +2 /
    /// mismatch −4 (A>G and T>C are charged like any mismatch and counted
    /// separately for the density gates), affine gaps `−(open + ext·n)`,
    /// introns by splice signal plus known-junction bonus minus the intron
    /// length penalty, soft clips free. N on either side scores 0 and counts
    /// as neither match nor mismatch. Identical constants to the extender
    /// and splicer, so scores are comparable across every producing path.
    fn rescore_strict(&mut self, query: &[u8], p: &mut Placement) {
        self.rescore_with(query, p, EditDir::Off);
    }

    /// [`rescore_strict`] with one editing direction exempt: substitutions of
    /// that direction score 0 instead of `mismatch` (they are still counted
    /// in `ea_ag`/`ea_tc`). This is the reduced-alphabet scoring of the
    /// hyperediting literature (Porath, Carmi & Levanon, Nat Commun 5:4726,
    /// 2014) applied to placement arbitration: a read whose strict-domain
    /// winner carries a hyperediting signature is re-arbitrated in the
    /// alphabet its own edits define, so that its true locus is not out-
    /// scored by a paralog that merely lacks the edits, and its edited tail
    /// is not soft-clipped away. The four-letter `unexplained` density cap in
    /// `mapq_of` still applies, so a paralog exploiting both directions or
    /// ordinary mismatches is not rescued by this.
    fn rescore_with(&mut self, query: &[u8], p: &mut Placement, dir: EditDir) {
        let ep = self.config.extend;
        let sp = SpliceParams::default();
        let contig = &self.index.reference.contigs[p.contig as usize];
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        let mut s = StrictScore {
            ..Default::default()
        };
        let mut score: i32 = 0;
        let mut ro = 0u32;
        let mut rf = p.pos;
        let mut scratch = std::mem::take(&mut self.scratch.rescore);
        let base = |i: u32| -> u8 {
            if i < contig.len {
                contig.base(i).to_ascii()
            } else {
                b'N'
            }
        };
        for op in &p.cigar {
            match op {
                CigarOp::Match(n) => {
                    let avail = (*n).min(contig.len.saturating_sub(rf));
                    scratch.clear();
                    if avail > 0 {
                        contig.decode_append(rf, rf + avail, &mut scratch);
                    }
                    for k in 0..*n as usize {
                        let qb = query.get(ro as usize + k).copied().unwrap_or(b'N');
                        let rb = if k < scratch.len() { scratch[k] } else { b'N' };
                        if qb == b'N' || rb == b'N' {
                            continue; // 0 score, uncounted
                        }
                        s.aligned += 1;
                        if qb == rb {
                            score += ep.match_score;
                            s.matches += 1;
                        } else {
                            if !dir.exempt(rb, qb) {
                                score += ep.mismatch;
                            }
                            match crate::extend::exempt_dir(rb, qb) {
                                Some(true) => s.ea_ag += 1,
                                Some(false) => s.ea_tc += 1,
                                None => s.mm += 1,
                            }
                        }
                    }
                    ro += *n;
                    rf += *n;
                }
                CigarOp::Ins(n) => {
                    score -= ep.gap_open + ep.gap_ext * *n as i32;
                    s.aligned += *n;
                    s.ins += 1;
                    ro += *n;
                }
                CigarOp::Del(n) => {
                    score -= ep.gap_open + ep.gap_ext * *n as i32;
                    s.aligned += *n;
                    s.del += 1;
                    rf += *n;
                }
                CigarOp::SoftClip(n) => {
                    s.clipped += *n;
                    ro += *n;
                }
                CigarOp::RefSkip(n) => {
                    let before = score;
                    // Annotated introns receive the known-junction bonus on
                    // either read orientation. Signal and length penalties
                    // apply to novel junctions only.
                    if lib.contains_any_strand(p.contig, rf, rf + *n) {
                        score += sp.known_bonus;
                    } else {
                        // Donor/acceptor are the first and last two bases
                        // OF the intron (GT..AG), matching `splice::dinuc`'s
                        // convention. The transcript's strand is unknown
                        // here — an unstranded library sequences half of
                        // its reads antisense — so a canonical motif on
                        // either strand counts and the read's own strand
                        // is not consulted.
                        let donor = [base(rf), base(rf.saturating_add(1))];
                        let acc = [
                            base(rf + n.saturating_sub(2)),
                            base(rf + n.saturating_sub(1)),
                        ];
                        let signal = splice_signal(donor, acc, false)
                            .score(&sp)
                            .max(splice_signal(donor, acc, true).score(&sp));
                        score += signal - intron_length_penalty(*n);
                    }
                    s.n_splice += 1;
                    s.splice_score += score - before;
                    rf += *n;
                }
            }
        }
        s.score = score;
        p.strict = s;
        p.ea_dir = dir;
        self.scratch.rescore = scratch;
    }

    /// Re-extend a soft-clipped genomic candidate with one editing direction
    /// exempt, over the window its current footprint plus clips define.
    /// Replaces `pos`/`cigar` when the new extension aligns more of the read.
    /// Spliced candidates are left alone (the window would have to span the
    /// intron; their clips are refined elsewhere).
    fn reextend_with(&mut self, query: &[u8], p: &mut Placement, dir: EditDir) {
        if p.cigar.iter().any(|op| matches!(op, CigarOp::RefSkip(_))) {
            return;
        }
        let read_len = query.len() as u32;
        let (mut lead, mut tail, mut ref_span, mut aligned) = (0u32, 0u32, 0u32, 0u32);
        let mut seen_m = false;
        for op in &p.cigar {
            match op {
                CigarOp::SoftClip(n) => {
                    if seen_m {
                        tail += *n;
                    } else {
                        lead += *n;
                    }
                }
                CigarOp::Match(n) => {
                    seen_m = true;
                    ref_span += *n;
                    aligned += *n;
                }
                CigarOp::Ins(n) => aligned += *n,
                CigarOp::Del(n) => ref_span += *n,
                CigarOp::RefSkip(_) => {}
            }
        }
        let contig = &self.index.reference.contigs[p.contig as usize];
        let wstart = p.pos.saturating_sub(lead + 30);
        let wend = (p.pos + ref_span + tail + 30).min(contig.len);
        if wend <= wstart {
            return;
        }
        let window = self.take_win(wstart, wend, contig);
        let hint = DiagHint {
            offset: (p.pos - wstart) as i64 - i64::from(lead),
            num: 1,
            den: 1,
        };
        let mut ep = self.config.extend;
        ep.edit_dir = dir;
        let ext = extend_hint(query, &window, &ep, &mut self.buf, hint);
        self.scratch.win = window;
        let new_aligned: u32 = ext
            .cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) => *n,
                _ => 0,
            })
            .sum();
        if ext.score > 0 && new_aligned >= aligned && new_aligned <= read_len {
            p.pos = wstart + ext.ref_start;
            p.cigar = ext.cigar;
            p.ext_score = ext.score;
        }
    }

    /// Count non-overlapping exact k-mers along the placement's M blocks
    /// (k from the seed config). Feeds `pen_cm` in the MAPQ formula for L1
    /// and genomic candidates alike — an anchor count measured on the
    /// placement, not borrowed from a chain.
    fn exact_kmer_support(&mut self, query: &[u8], p: &Placement) -> usize {
        let k = self.config.seed.k as usize;
        if k == 0 || query.is_empty() {
            return 0;
        }
        let contig = &self.index.reference.contigs[p.contig as usize];
        let mut count = 0usize;
        let mut ro = 0usize;
        let mut rf = p.pos;
        let mut scratch = std::mem::take(&mut self.scratch.kmer_scan);
        for op in &p.cigar {
            match op {
                CigarOp::Match(n) => {
                    let n = *n as usize;
                    let avail = (n as u32).min(contig.len.saturating_sub(rf));
                    scratch.clear();
                    if avail > 0 {
                        contig.decode_append(rf, rf + avail, &mut scratch);
                    }
                    let mut i = 0usize;
                    while i + k <= n {
                        let mut ok = true;
                        for t in 0..k {
                            let qb = query.get(ro + i + t).copied().unwrap_or(b'N');
                            let rb = if i + t < scratch.len() {
                                scratch[i + t]
                            } else {
                                b'N'
                            };
                            // A k-mer is exact in the alphabet the candidate
                            // was scored in: for a hyperedited read that is
                            // its single-direction reduced alphabet, where its
                            // edits are identities. Counting four-letter exact
                            // k-mers on such a read gives structurally ~0
                            // anchors (30 edits in 150 bp break every 15-mer)
                            // and hence MAPQ 0 whatever the margin.
                            if qb != rb && !p.ea_dir.exempt(rb, qb) {
                                ok = false;
                                break;
                            }
                        }
                        if ok {
                            count += 1;
                            i += k;
                        } else {
                            i += 1;
                        }
                    }
                    ro += n;
                    rf += n as u32;
                }
                CigarOp::Ins(n) | CigarOp::SoftClip(n) => ro += *n as usize,
                CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
            }
        }
        self.scratch.kmer_scan = scratch;
        count
    }

    /// Genomic layer driver, in candidate form: seed → anchor → chain →
    /// build a locus-deduplicated chain list → extend each chain
    /// (intron-chain / spliced / EA-Myers / fast lane / banded DP) into
    /// 0..2 placements → split/span rescue on the best placement. Each
    /// defensible locus enters arbitration in [`Aligner::align_read`].
    fn genomic_candidates(&mut self, seq: &[u8]) -> GenomicBatch {
        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;

        let mut mins = std::mem::take(&mut self.scratch.mins_out);
        minimizers_into(seq, cfg.seed, &mut self.scratch.seed, &mut mins);
        if mins.is_empty() {
            self.scratch.mins_out = mins;
            return GenomicBatch {
                placements: Vec::new(),
                n_seeds: 0,
                seed_budget: None,
            };
        }

        // Round 1: default-density seeding, occurrence cap = freq_cutoff.
        let r1 = self.chain_round(&mins, read_len, k, self.index.freq_cutoff);

        // Round 2 raises the occurrence cap when the first round covers
        // little of the read, retaining candidates from high-copy repeats.
        // The anchor limit bounds chaining cost on repetitive sequence.
        let r1_cov = r1
            .as_ref()
            .and_then(|r| r.best())
            .map(|b| {
                let (qs, qe) = b.read_span(k);
                f64::from(qe.saturating_sub(qs)) / f64::from(read_len.max(1))
            })
            .unwrap_or(0.0);
        let mut seed_budget = None;
        let mut r2 = if r1_cov < MID_ROUND_MIN_COVERAGE {
            let (round, budget) = self.bounded_chain_round(&mins, read_len, k);
            if budget.limited() {
                seed_budget = Some(budget);
            }
            round
        } else {
            None
        };

        // Truncating technical tails can change terminal minimizer selection
        // even when the necessary k-mer lies wholly inside the retained query.
        // Only a complete chaining failure permits query-edge completion;
        // successful rounds, DNA and edit-variant search are left untouched.
        let ordinary_failed = cfg.rna
            && !cfg.extend.edit_dir.is_on()
            && r1.as_ref().is_none_or(|r| r.best().is_none())
            && r2.as_ref().is_none_or(|r| r.best().is_none());
        if ordinary_failed
            && append_query_edge_kmers(seq.len(), cfg.seed, &self.scratch.seed, &mut mins) > 0
        {
            if let Some(failed) = r2.take() {
                self.scratch.chains_out = failed.chains;
            }
            let (round, budget) = self.bounded_chain_round(&mins, read_len, k);
            if budget.limited() {
                seed_budget = Some(budget);
            }
            r2 = round;
        }

        // Editing-aware variant round: only when chaining fully failed (weak).
        let weak = r1
            .as_ref()
            .and_then(|r| r.best())
            .map(|b| b.score)
            .unwrap_or(0)
            == 0;
        let r3 = if cfg.extend.edit_dir.is_on() && weak {
            // Base anchors first: edit-variant anchors must not alter a clean chain.
            let mut base_anchors = std::mem::take(&mut self.scratch.anchors);
            collect_anchors_into(
                self.index,
                &mins,
                read_len,
                k,
                MID_OCC_CAP,
                &mut base_anchors,
            );
            if !base_anchors.is_empty() {
                let mut cp = cfg.chain;
                cp.min_chain_score = 25;
                let mut chains = std::mem::take(&mut self.scratch.chains_out);
                chain_anchors_into(&mut base_anchors, &cp, &mut chains, &mut self.scratch.chain);
                self.scratch.anchors = base_anchors;
                let second = second_score_with(&chains, &mut self.scratch.second_keys);
                Some(RoundOutcome { chains, second })
            } else {
                self.scratch.anchors = base_anchors;
                let (mut anchors, _var_hits) = crate::index::collect_anchors_edit_variants(
                    self.index,
                    &mins,
                    read_len,
                    k,
                    MID_OCC_CAP,
                );
                if anchors.is_empty() {
                    None
                } else {
                    let mut cp = cfg.chain;
                    cp.min_chain_score = 25;
                    let mut chains = std::mem::take(&mut self.scratch.chains_out);
                    chain_anchors_into(&mut anchors, &cp, &mut chains, &mut self.scratch.chain);
                    let second = second_score_with(&chains, &mut self.scratch.second_keys);
                    Some(RoundOutcome { chains, second })
                }
            }
        } else {
            None
        };

        let mut placements = self.extend_genomic_rounds(
            seq,
            &mins,
            r1.as_ref(),
            r2.as_ref(),
            r3.as_ref(),
            seed_budget.is_some(),
        );
        if let Some(r1) = r1 {
            self.scratch.chains_out = r1.chains;
        }
        // Split/span rescue and micro-deletion conversion on the best
        // genomic placement (by extension score, chain order on ties) — the
        // same post-processing the single-winner path applied to its winner.
        let mut best_idx: Option<usize> = None;
        for (i, p) in placements.iter().enumerate() {
            let better = match best_idx {
                None => true,
                Some(j) => {
                    p.ext_score > placements[j].ext_score
                        || (p.ext_score == placements[j].ext_score
                            && p.chain_score > placements[j].chain_score)
                }
            };
            if better {
                best_idx = Some(i);
            }
        }
        if let Some(i) = best_idx {
            let p = placements[i].clone();
            let tmp = ReadAlignment {
                contig: p.contig,
                pos: p.pos,
                strand: p.strand,
                score: p.ext_score,
                chain_score: p.chain_score,
                cigar: p.cigar,
                junctions: p.junctions,
                ..Default::default()
            };
            let tmp = self.apply_split_rescue(seq, tmp);
            let tmp = self.apply_span_rescue(seq, tmp);
            let tmp = self.convert_annotated_micro_deletions(tmp);
            let q = &mut placements[i];
            q.pos = tmp.pos;
            q.cigar = tmp.cigar;
            q.junctions = tmp.junctions;
        }

        // Additional context candidates must not replace or consume the old
        // chain/splice/indel extension work above. A perfect whole-query edge
        // hit attains the ordinary score ceiling and bypasses this search.
        if cfg.repeat_context
            && ordinary_failed
            && !placements.iter().any(|p| {
                p.ext_score == cfg.extend.match_score * read_len as i32
                    && p.cigar == [CigarOp::Match(read_len)]
            })
        {
            let context = if cfg.context_background {
                crate::repeat_context::generate(self.index, seq, cfg.extend)
            } else {
                crate::repeat_context::generate_with_background(self.index, seq, cfg.extend, false)
            };
            if let Some((proposals, context)) = context {
                seed_budget
                    .get_or_insert_with(Default::default)
                    .repeat_context = Some(context);
                placements.extend(proposals.into_iter().map(|p| {
                    let mut candidate = Placement::from_genomic(
                        p.contig,
                        p.pos,
                        p.strand,
                        p.cigar,
                        Vec::new(),
                        0,
                        p.score,
                    );
                    candidate.origin = crate::mapq::Origin::RepeatContext;
                    candidate
                }));
            }
        }
        let n_seeds = mins.len();
        self.scratch.mins_out = mins;

        GenomicBatch {
            placements,
            n_seeds,
            seed_budget,
        }
    }

    /// Schedule and extend genomic chains with one shared DP budget.
    fn extend_genomic_rounds(
        &mut self,
        seq: &[u8],
        mins: &[Minimizer],
        r1: Option<&RoundOutcome>,
        r2: Option<&RoundOutcome>,
        r3: Option<&RoundOutcome>,
        seed_search_limited: bool,
    ) -> Vec<Placement> {
        let k = self.config.seed.k;
        // A limited repeat search may introduce higher-chain-score competitors
        // that exhaust extension before a complete first-round path is tried.
        // Preserve that round's existing budget priority; additional proposals
        // still compete by common sequence score. Complete search is unchanged.
        let rounds: Vec<&RoundOutcome> = r1.into_iter().chain(r2).chain(r3).collect();
        let mut win_idx: Option<usize> = None;
        for (i, r) in rounds.iter().enumerate() {
            let better = match (win_idx, r.best()) {
                (None, Some(_)) => true,
                (Some(w), Some(b)) => b.score > rounds[w].best().map(|x| x.score).unwrap_or(0),
                _ => false,
            };
            if better {
                win_idx = Some(i);
            }
        }
        if seed_search_limited && r1.is_some_and(|r| r.best().is_some()) {
            win_idx = Some(0);
        }
        let Some(w) = win_idx else {
            return Vec::new();
        };
        let cr = rounds[w];
        let second = cr.second;
        let extra_top: Vec<Chain> = rounds
            .iter()
            .flat_map(|r| r.chains.iter().take(TOP_CANDIDATES).cloned())
            .collect();
        let best = match cr.best() {
            Some(b) => b.clone(),
            None => {
                return Vec::new();
            }
        };

        // Candidate chains: deduplicated by locus (same contig+strand,
        // ref-span start within 50bp), in chain-score order.
        let mut cand_chains: Vec<Chain> = Vec::new();
        let push_dedup = |c: Chain, out: &mut Vec<Chain>| {
            let dup = out.iter().any(|x| {
                x.contig == c.contig
                    && x.strand == c.strand
                    && x.ref_span(k).0.abs_diff(c.ref_span(k).0) < 50
            });
            if !dup {
                out.push(c);
            }
        };
        for c in cr.chains.iter().take(TOP_CANDIDATES) {
            push_dedup(c.clone(), &mut cand_chains);
        }
        for c in extra_top {
            if cand_chains.len() >= MAX_TIED_CANDIDATES {
                break;
            }
            push_dedup(c, &mut cand_chains);
        }

        let tied = tied_best(&cr.chains);
        let decisive = tied.len() <= 1 && second * 5 < best.score * 4;
        let rc_storage = self.take_rc(seq);
        let mut placements: Vec<Placement> = Vec::new();
        let mut dp_budget = EXTEND_BUDGET;
        for (i, chain) in cand_chains.iter().take(MAX_TIED_CANDIDATES).enumerate() {
            let is_best = i == 0;
            self.extend_chain(
                seq,
                mins,
                chain,
                is_best,
                decisive && is_best,
                &rc_storage,
                &mut dp_budget,
                seed_search_limited,
                &mut placements,
            );
        }
        self.scratch.rc = rc_storage;
        placements
    }

    /// Extend one chain into spliced and continuous placement candidates.
    /// RNA chains also retain segment-wise continuous explanations. The
    /// unspliced extension uses EA-Myers with banded DP while budget lasts.
    /// `fast_lane_ok` gates the Kadane fast lane (decisive best chain only).
    #[allow(clippy::too_many_arguments)]
    fn extend_chain(
        &mut self,
        seq: &[u8],
        mins: &[Minimizer],
        chain: &Chain,
        is_best: bool,
        fast_lane_ok: bool,
        rc_storage: &[u8],
        dp_budget: &mut usize,
        seed_search_limited: bool,
        out: &mut Vec<Placement>,
    ) {
        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;

        // Spliced interpretation (RNA; splice gate: <2 segments
        // has no spliced reading).
        if cfg.rna {
            let min_intron = SpliceParams::default().min_intron;
            let segs = segment_chain(chain, k, min_intron);
            if segs.len() >= 2 {
                // A repeat-derived anchor can attach to a genuine continuous
                // placement and manufacture a plausible intron. A successful
                // splice must not suppress the continuous explanation at any
                // of its collinear segments. Longer segments get the shared
                // extension budget first.
                let mut ordered = segs.clone();
                ordered.sort_by_key(|s| std::cmp::Reverse(s.q_end - s.q_start));
                for seg in ordered {
                    let segment = Chain {
                        contig: chain.contig,
                        strand: chain.strand,
                        score: chain.score,
                        anchors: chain
                            .anchors
                            .iter()
                            .copied()
                            .filter(|a| {
                                a.qpos >= seg.q_start
                                    && a.qpos + k <= seg.q_end
                                    && a.rpos >= seg.r_start
                                    && a.rpos + k <= seg.r_end
                            })
                            .collect(),
                    };
                    if !segment.anchors.is_empty() {
                        self.extend_unspliced_chain(
                            seq, &segment, false, rc_storage, dp_budget, out,
                        );
                    }
                }
            }
            if segs.len() >= 2 && is_best {
                // The intron fast path may choose any exon pair supported by
                // the complete anchor pool; extra segments in this chain do
                // not invalidate that pair. Keep its own boundary geometry.
                let empty_lib = JunctionLib::default();
                let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
                let query: &[u8] = match chain.strand {
                    Strand::Plus => seq,
                    // Same bytes as a fresh revcomp(seq): the caller passes
                    // exactly that buffer.
                    Strand::Minus => rc_storage,
                };
                let mut ic_anchors = std::mem::take(&mut self.scratch.anchors);
                if seed_search_limited {
                    crate::index::collect_anchors_budgeted_into(
                        self.index,
                        mins,
                        read_len,
                        k,
                        MID_OCC_CAP,
                        MID_ROUND_MAX_ANCHORS,
                        &mut ic_anchors,
                    );
                } else {
                    collect_anchors_into(
                        self.index,
                        mins,
                        read_len,
                        k,
                        MID_OCC_CAP,
                        &mut ic_anchors,
                    );
                }
                ic_anchors.retain(|a| a.contig == chain.contig && a.strand == chain.strand);
                let ic_hit = crate::intron_chain::try_intron_chain_placement(
                    self.index.reference,
                    lib,
                    query,
                    &ic_anchors,
                    k,
                    &crate::intron_chain::IntronParams::default(),
                    &cfg.extend,
                    &mut self.buf,
                );
                self.scratch.anchors = ic_anchors;
                if let Some(hit) = ic_hit {
                    let placement = Placement::from_genomic(
                        chain.contig,
                        hit.genomic_start(),
                        chain.strand,
                        hit.cigar,
                        hit.junctions,
                        chain.score,
                        hit.extension.score,
                    );
                    // The fast path searches the whole anchor set and may
                    // choose a different exon pair from this chain. Its CIGAR
                    // cannot be refined by rewriting N with this chain's
                    // junction. Keep only internally consistent proposals.
                    if placement.has_consistent_junctions() {
                        out.push(placement);
                    }
                    // A fast proposal never suppresses the generic spliced
                    // explanation; common rescoring decides between them.
                }
            }
            if segs.len() >= 2 {
                if let Some(hit) = self.try_spliced(seq, rc_storage, chain) {
                    out.push(Placement::from_genomic(
                        chain.contig,
                        hit.pos,
                        chain.strand,
                        hit.cigar,
                        hit.junctions,
                        chain.score,
                        hit.score,
                    ));
                    return;
                }
            }
        }

        self.extend_unspliced_chain(seq, chain, fast_lane_ok, rc_storage, dp_budget, out);
    }

    /// Continuous explanation of one collinear chain, without splice search.
    fn extend_unspliced_chain(
        &mut self,
        seq: &[u8],
        chain: &Chain,
        fast_lane_ok: bool,
        rc_storage: &[u8],
        dp_budget: &mut usize,
        out: &mut Vec<Placement>,
    ) {
        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;
        // Unspliced interpretation — fast lane (decisive best chain only:
        // zero-mismatch Kadane, point-equivalent to the affine DP).
        if fast_lane_ok {
            if let Some((ext, wstart)) = self.try_fast_lane(seq, rc_storage, chain, read_len, k) {
                if ext.score > 0 {
                    out.push(Placement::from_genomic(
                        chain.contig,
                        wstart + ext.ref_start,
                        chain.strand,
                        ext.cigar,
                        Vec::new(),
                        chain.score,
                        ext.score,
                    ));
                    return;
                }
            }
        }

        // Unspliced interpretation — EA-Myers window extend (≤256bp; the
        // bit-parallel infix search with Kadane-optimal soft clips).
        let mut ea_done_here = false;
        if read_len > 0 && read_len <= 256 {
            let dir = cfg.extend.edit_dir;
            let m = read_len as usize;
            let match_s = cfg.extend.match_score;
            let mm_s = cfg.extend.mismatch;
            let threshold = ((read_len as i32) / 3).max(10);
            let contig = &self.index.reference.contigs[chain.contig as usize];
            let (wstart, wend) = continuous_extension_window(chain, k, read_len, contig.len);
            let window = self.take_win(wstart, wend, contig);
            let q: &[u8] = match chain.strand {
                Strand::Plus => seq,
                Strand::Minus => rc_storage,
            };
            if window.len() >= m {
                let (mut myers_a, mut myers_b) = (
                    std::mem::take(&mut self.scratch.myers_a),
                    std::mem::take(&mut self.scratch.myers_b),
                );
                let (d, start) = if m <= 128 {
                    crate::myers_ea::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
                } else {
                    crate::myers_ea::long::infix_best_start_into(
                        q,
                        &window,
                        &mut myers_a,
                        &mut myers_b,
                    )
                };
                self.scratch.myers_a = myers_a;
                self.scratch.myers_b = myers_b;
                let mut usable = d <= threshold && start + m <= window.len();
                if usable {
                    // Per-base verification: mm must exactly explain d,
                    // otherwise an indel is indicated.
                    let mut mm_full: i32 = 0;
                    let mut has_n = false;
                    for kk in 0..m {
                        let qb = q[kk];
                        let rb = window[start + kk];
                        if qb == b'N' || rb == b'N' {
                            has_n = true;
                            break;
                        }
                        if qb != rb && crate::extend::exempt_dir(rb, qb).is_none() {
                            mm_full += 1;
                        }
                    }
                    if has_n {
                        usable = false;
                    } else if mm_full != d {
                        usable = false; // indel implied → DP below
                    }
                }
                if usable {
                    // Optimal soft clip ≡ Kadane max subarray on the diagonal.
                    let mut best_sum = 0i32;
                    let (mut bi, mut bj) = (0usize, 0usize);
                    let mut cur_sum = 0i32;
                    let mut cur_i = 0usize;
                    for kk in 0..m {
                        let qb = q[kk];
                        let rb = window[start + kk];
                        let s = if qb == rb {
                            match_s
                        } else if dir.exempt(rb, qb) {
                            0
                        } else {
                            mm_s
                        };
                        if cur_sum <= 0 {
                            cur_sum = s;
                            cur_i = kk;
                        } else {
                            cur_sum += s;
                        }
                        if cur_sum > best_sum {
                            best_sum = cur_sum;
                            bi = cur_i;
                            bj = kk + 1;
                        }
                    }
                    let seg_len = (bj - bi) as u32;
                    if best_sum > 0 && seg_len * 2 >= read_len {
                        let mut cigar = Vec::with_capacity(3);
                        if bi > 0 {
                            cigar.push(CigarOp::SoftClip(bi as u32));
                        }
                        cigar.push(CigarOp::Match(seg_len));
                        if m - bj > 0 {
                            cigar.push(CigarOp::SoftClip((m - bj) as u32));
                        }
                        let pos = wstart + (start + bi) as u32;
                        out.push(Placement::from_genomic(
                            chain.contig,
                            pos,
                            chain.strand,
                            cigar,
                            Vec::new(),
                            chain.score,
                            best_sum,
                        ));
                        ea_done_here = true;
                    }
                }
            }
            self.scratch.win = window;
        }
        if ea_done_here {
            return;
        }
        if ea_done_here {
            return;
        }

        // Unspliced interpretation — banded DP (indel fallback and long
        // reads), within the shared extension budget.
        if *dp_budget == 0 {
            return;
        }
        *dp_budget -= 1;
        let contig = &self.index.reference.contigs[chain.contig as usize];
        let (rs, _) = chain.ref_span(k);
        let (wstart, wend) = continuous_extension_window(chain, k, read_len, contig.len);
        let window = self.take_win(wstart, wend, contig);
        let q: &[u8] = match chain.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };
        let (qs, _) = chain.read_span(k);
        let hint = DiagHint {
            offset: (rs - wstart) as i64 - qs as i64,
            num: 1,
            den: 1,
        };
        let ext = extend_hint(q, &window, &cfg.extend, &mut self.buf, hint);
        self.scratch.win = window;
        if ext.score > 0 {
            out.push(Placement::from_genomic(
                chain.contig,
                wstart + ext.ref_start,
                chain.strand,
                ext.cigar,
                Vec::new(),
                chain.score,
                ext.score,
            ));
        }
    }

    /// Splice path: segment → junction refine → pseudo-reference extend →
    /// N CIGAR. splice gate: <2 segments ⇒ None (align_spliced would
    /// return None anyway — skip the revcomp + entry cost).
    fn try_spliced(
        &mut self,
        seq: &[u8],
        rc_storage: &[u8],
        best: &Chain,
    ) -> Option<ReadAlignment> {
        let min_intron = SpliceParams::default().min_intron;
        let segs = segment_chain(best, self.config.seed.k, min_intron);
        if segs.len() < 2 {
            return None;
        }
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        // rc_storage is the caller's revcomp(seq) buffer (same bytes).
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };
        let hit = align_spliced(
            self.index.reference,
            lib,
            query,
            best,
            self.config.seed.k,
            &SpliceParams::default(),
            &self.config.extend,
            &mut self.buf,
        )?;
        if hit.extension.score <= 0 {
            return None;
        }
        let mut aln = ReadAlignment {
            contig: best.contig,
            pos: hit.genomic_start(),
            strand: best.strand,
            score: hit.extension.score,
            chain_score: best.score,
            second_chain_score: 0,
            cigar: hit.cigar,
            n_anchors: best.anchors.len(),
            junctions: hit.junctions,
            ..Default::default()
        };
        self.recount_mm_ea(query, &mut aln);
        Some(aln)
    }

    /// Compare rescue proposals with the primary in the same strict domain.
    /// Extension totals and genomic-scan ranks are not read-placement scores.
    /// The existing split cost is an acceptance margin; a library bonus alone
    /// must not rewrite an equally good continuous sequence explanation.
    fn select_rescue(
        &mut self,
        seq: &[u8],
        original: ReadAlignment,
        mut proposal: ReadAlignment,
    ) -> ReadAlignment {
        if original.pos == proposal.pos && original.cigar == proposal.cigar {
            return original;
        }
        let mut before = Placement::from_genomic(
            original.contig,
            original.pos,
            original.strand,
            original.cigar.clone(),
            original.junctions.clone(),
            original.chain_score,
            original.score,
        );
        let mut after = Placement::from_genomic(
            proposal.contig,
            proposal.pos,
            proposal.strand,
            proposal.cigar.clone(),
            proposal.junctions.clone(),
            proposal.chain_score,
            proposal.score,
        );
        let rc_storage = self.take_rc(seq);
        let query = if original.strand == Strand::Minus {
            &rc_storage
        } else {
            seq
        };
        self.rescore_strict(query, &mut before);
        self.rescore_strict(query, &mut after);
        self.scratch.rc = rc_storage;
        let sequence_gain = after.strict.score - after.strict.splice_score
            > before.strict.score - before.strict.splice_score;
        if sequence_gain && after.strict.score > before.strict.score + crate::splice::SPLIT_COST {
            proposal.score = after.strict.score;
            proposal
        } else {
            original
        }
    }

    /// Split rescue (both tails), accepted against the original primary.
    fn apply_split_rescue(&mut self, seq: &[u8], aln: ReadAlignment) -> ReadAlignment {
        let proposal = self.propose_split_rescue(seq, aln.clone());
        self.select_rescue(seq, aln, proposal)
    }

    fn propose_split_rescue(&mut self, seq: &[u8], aln: ReadAlignment) -> ReadAlignment {
        let empty_lib = JunctionLib::default();
        // Query in extension orientation: (RC'd if minus) read; tails operate
        // in extension coordinates; result CIGAR is plus-reference isomorphic.
        let rc_storage = self.take_rc(seq);
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let mut rs = 0u32;
        let mut re = query.len() as u32;
        let mut seen_body = false;
        for op in &aln.cigar {
            match op {
                CigarOp::SoftClip(n) => {
                    if !seen_body {
                        rs = *n;
                    } else {
                        re -= *n;
                    }
                }
                _ => seen_body = true,
            }
        }
        let read_start = rs;
        let read_end = re;
        let mut out = aln;
        let primary_ref_end = {
            let mut e = out.pos;
            for op in &out.cigar {
                match op {
                    CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => e += *n,
                    _ => {}
                }
            }
            e
        };
        let params = SpliceParams::default();
        // Right tail = high-coordinate side (downstream intron), left tail =
        // low side; same on both strands in extension coordinates.
        if query.len() as u32 - read_end >= 1 {
            let r_tail = query.len() as u32 - read_end;
            let ctx = crate::split::SplitContext {
                reference: self.index.reference,
                index: Some(self.index),
                lib,
                read: query,
                contig: out.contig,
                strand: out.strand,
                pos: out.pos,
                ref_end: primary_ref_end,
                read_start,
                read_end,
                cigar: &out.cigar,
                seed_params: self.config.seed,
                splice_params: params,
                extend_params: self.config.extend,
            };
            if let Some(rescue) = crate::split::rescue_right_tail(&ctx, &mut self.buf) {
                if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                    // tail-rescue gate: tails <12bp require either
                    // known_support ≥ 2 OR (canonical signal + zero tail-vs-ref
                    // mismatch) — a de-novo short tail with a perfect anchor-side
                    // junction match is accepted.
                    let t3_reject = r_tail < 12 && {
                        let contig = &self.index.reference.contigs[out.contig as usize];
                        let mm =
                            count_tail_mismatches(&rescue.cigar, rescue.pos, contig, query, true);
                        if rescue.junction.known_support < 2 {
                            mm > 0 || rescue.junction.signal == SpliceSignal::NonCanonical
                        } else {
                            mm > 0
                        }
                    };
                    if !t3_reject {
                        out.pos = rescue.pos;
                        out.cigar = rescue.cigar;
                        out.junctions.push(rescue.junction);
                        out.score += rescue.score;
                        self.scratch.rc = rc_storage;
                        return out;
                    }
                }
            }
        }
        if read_start >= 1 {
            let ctx = crate::split::SplitContext {
                reference: self.index.reference,
                index: Some(self.index),
                lib,
                read: query,
                contig: out.contig,
                strand: out.strand,
                pos: out.pos,
                ref_end: primary_ref_end,
                read_start,
                read_end,
                cigar: &out.cigar,
                seed_params: self.config.seed,
                splice_params: params,
                extend_params: self.config.extend,
            };
            if let Some(rescue) = crate::split::rescue_left_tail(&ctx, &mut self.buf) {
                if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                    // Mirror of the right-tail gate:
                    // de-novo short tails accepted with canonical signal +
                    // zero tail-vs-ref mismatch.
                    let t3_reject = read_start < 12 && {
                        let contig = &self.index.reference.contigs[out.contig as usize];
                        let mm =
                            count_tail_mismatches(&rescue.cigar, rescue.pos, contig, query, false);
                        if rescue.junction.known_support < 2 {
                            mm > 0 || rescue.junction.signal == SpliceSignal::NonCanonical
                        } else {
                            mm > 0
                        }
                    };
                    if !t3_reject {
                        out.pos = rescue.pos;
                        out.cigar = rescue.cigar;
                        out.junctions.push(rescue.junction);
                        out.score += rescue.score;
                    }
                }
            }
        }
        self.scratch.rc = rc_storage;
        out
    }

    /// Library-driven run-through reinterpretation (span.rs; skipped without
    /// a library or when the alignment is already spliced).
    fn apply_span_rescue(&mut self, seq: &[u8], aln: ReadAlignment) -> ReadAlignment {
        let proposal = self.propose_span_rescue(seq, aln.clone());
        self.select_rescue(seq, aln, proposal)
    }

    fn propose_span_rescue(&mut self, seq: &[u8], mut aln: ReadAlignment) -> ReadAlignment {
        let rc_storage = self.take_rc(seq);
        let Some(lib) = self.jlib.as_deref() else {
            self.scratch.rc = rc_storage;
            return aln;
        };
        if !aln.junctions.is_empty() {
            self.scratch.rc = rc_storage;
            return aln; // already has a spliced interpretation
        };
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let ctx = crate::split::SplitContext {
            reference: self.index.reference,
            index: Some(self.index),
            lib,
            read: query,
            contig: aln.contig,
            strand: aln.strand,
            pos: aln.pos,
            ref_end: 0, // span recomputes ref_end from pos+cigar
            read_start: 0,
            read_end: 0,
            cigar: &aln.cigar,
            seed_params: self.config.seed,
            splice_params: SpliceParams::default(),
            extend_params: self.config.extend,
        };
        if let Some(rescue) = crate::span::rescue_span(&ctx, &mut self.buf) {
            if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                aln.cigar = rescue.cigar;
                aln.junctions.push(rescue.junction);
                aln.score += rescue.score;
            }
        }
        self.scratch.rc = rc_storage;
        aln
    }

    /// Annotated micro-deletion → junction: a D of ≤50bp exactly matching a
    /// library junction is rewritten to N (1-40bp annotated introns and
    /// deletions are sequence-indistinguishable; only the library tells).
    fn convert_annotated_micro_deletions(&self, mut aln: ReadAlignment) -> ReadAlignment {
        let Some(lib) = self.jlib.as_deref() else {
            return aln;
        };
        let mut pos = aln.pos;
        let mut new_cigar = Vec::with_capacity(aln.cigar.len());
        for op in &aln.cigar {
            let adv = match op {
                CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => *n,
                _ => 0,
            };
            let converted = match op {
                CigarOp::Del(n) if *n <= 50 => {
                    let j = Junction {
                        contig: aln.contig,
                        start: pos,
                        end: pos + *n,
                        minus_strand: aln.strand == Strand::Minus,
                    };
                    if lib.contains(&j) {
                        aln.junctions.push(RefinedJunction {
                            junction: j,
                            signal: SpliceSignal::GtAg,
                            known_support: lib.support(&j),
                        });
                        CigarOp::RefSkip(*n)
                    } else {
                        *op
                    }
                }
                _ => *op,
            };
            new_cigar.push(converted);
            pos += adv;
        }
        aln.cigar = new_cigar;
        aln
    }

    /// No-DP fast lane for decisive, single-diagonal, fully-covered chains.
    ///
    /// Equivalence (zero-mismatch gate + Kadane max subarray): gap-free affine
    /// DP on one diagonal = Kadane; with zero mismatches the score is the
    /// theoretical maximum, so DP can improve neither via gaps nor off-diagonal.
    /// Returns (Extension, wstart) matching extend_hint's coordinates.
    fn try_fast_lane(
        &mut self,
        seq: &[u8],
        rc_storage: &[u8],
        best: &Chain,
        read_len: u32,
        k: u32,
    ) -> Option<(Extension, u32)> {
        let cfg = self.config;

        // Condition 2: all anchors on the same diagonal (ref_pos − read_pos).
        let first = best.anchors.first()?;
        let diag = first.rpos as i64 - first.qpos as i64;
        if !best
            .anchors
            .iter()
            .all(|a| a.rpos as i64 - a.qpos as i64 == diag)
        {
            return None;
        }

        // Condition 3: adjacent anchor read-gap ≤ k+w+5 (=25bp).
        let max_gap = k + cfg.seed.w + 5;
        for w in best.anchors.windows(2) {
            if w[1].qpos - w[0].qpos > max_gap {
                return None;
            }
        }

        // Condition 4: chain read span covers the full read (≤2bp overhangs).
        let (qstart, qend) = best.read_span(k);
        if qstart > 2 || read_len - qend > 2 {
            return None;
        }

        // rc_storage is the caller's revcomp(seq) buffer (same bytes).
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };

        // Window formula identical to the DP extend path.
        let contig = &self.index.reference.contigs[best.contig as usize];
        let (rs, _re) = best.ref_span(k);
        let flank = read_len / 2 + 30;
        let wstart = rs.saturating_sub(flank);
        let ep = cfg.extend;

        // Kadane scan along the diagonal. read_pos i ↔ contig_pos (i + diag).
        // Per-position scores match SubstLut: match=+2, mismatch=−4, EA=0, N=0.
        let mut best_score = 0i32;
        let mut best_start = 0usize;
        let mut best_end = 0usize;
        let mut cur_score = 0i32;
        let mut cur_start = 0usize;
        let mut cur_mm = 0u32;
        let mut best_mm = 0u32;

        for (i, &qb) in query.iter().enumerate() {
            let rpos = i as i64 + diag;
            let (s, is_mm) = if rpos < 0 || rpos >= contig.len as i64 {
                (0, false)
            } else {
                let rb = contig.base(rpos as u32).to_ascii();
                if qb == b'N' || rb == b'N' {
                    (0, false)
                } else if qb == rb {
                    (ep.match_score, false)
                } else if ep.edit_dir.exempt(rb, qb) {
                    (0, false)
                } else {
                    (ep.mismatch, true)
                }
            };
            // Kadane with ≤0 reset (matches DP H[i] = max(0, H[i−1] + s(i))).
            if cur_score <= 0 {
                cur_score = s;
                cur_start = i;
                cur_mm = if is_mm { 1 } else { 0 };
            } else {
                cur_score += s;
                if is_mm {
                    cur_mm += 1;
                }
            }
            if cur_score > best_score {
                best_score = cur_score;
                best_start = cur_start;
                best_end = i + 1;
                best_mm = cur_mm;
            }
        }

        if best_score <= 0 {
            return None;
        }

        // Condition 5 (zero-mismatch gate): any mismatch in the best subarray
        // opens a gap-improvement opportunity for DP → refuse the fast lane.
        if best_mm > 0 {
            return None;
        }

        let read_start = best_start as u32;
        let read_end = best_end as u32;
        let match_len = read_end - read_start;

        let mut cigar = Vec::with_capacity(3);
        if read_start > 0 {
            cigar.push(CigarOp::SoftClip(read_start));
        }
        cigar.push(CigarOp::Match(match_len));
        let trailing = read_len - read_end;
        if trailing > 0 {
            cigar.push(CigarOp::SoftClip(trailing));
        }

        let ref_start = (read_start as i64 + diag - wstart as i64) as u32;
        let ext = Extension {
            score: best_score,
            ref_start,
            read_start,
            read_end,
            cigar,
        };

        Some((ext, wstart))
    }

    /// mm/ea recount along the CIGAR, base-by-base (extension coordinates).
    /// ea = (ref A, read G)/(ref T, read C) editing-type mismatches; mm = rest.
    pub(crate) fn recount_mm_ea(&mut self, query: &[u8], aln: &mut ReadAlignment) {
        recount_mm_ea_in(self.index.reference, query, aln);
    }

    fn bounded_chain_round(
        &mut self,
        mins: &[Minimizer],
        read_len: u32,
        k: u32,
    ) -> (Option<RoundOutcome>, crate::index::AnchorBudget) {
        let mut anchors = std::mem::take(&mut self.scratch.anchors);
        let budget = crate::index::collect_anchors_budgeted_into(
            self.index,
            mins,
            read_len,
            k,
            MID_OCC_CAP,
            MID_ROUND_MAX_ANCHORS,
            &mut anchors,
        );
        if anchors.is_empty() {
            self.scratch.anchors = anchors;
            return (None, budget);
        }
        let mut chains = std::mem::take(&mut self.scratch.chains_out);
        chain_anchors_into(
            &mut anchors,
            &self.config.chain,
            &mut chains,
            &mut self.scratch.chain,
        );
        self.scratch.anchors = anchors;
        let second = second_score_with(&chains, &mut self.scratch.second_keys);
        (Some(RoundOutcome { chains, second }), budget)
    }

    pub(crate) fn chain_round(
        &mut self,
        mins: &[Minimizer],
        read_len: u32,
        k: u32,
        occ_cap: u32,
    ) -> Option<RoundOutcome> {
        let mut anchors = std::mem::take(&mut self.scratch.anchors);
        collect_anchors_into(self.index, mins, read_len, k, occ_cap, &mut anchors);
        if anchors.is_empty() {
            self.scratch.anchors = anchors;
            return None;
        }
        let mut chains = std::mem::take(&mut self.scratch.chains_out);
        chain_anchors_into(
            &mut anchors,
            &self.config.chain,
            &mut chains,
            &mut self.scratch.chain,
        );
        self.scratch.anchors = anchors;
        let second = second_score_with(&chains, &mut self.scratch.second_keys);
        Some(RoundOutcome { chains, second })
    }
}

/// Set `aln.mm_count` / `aln.ea_count` by walking the CIGAR base-by-base
/// against `reference`. `query` must already be in the alignment's frame
/// (reverse-complemented for [`Strand::Minus`]). Callers that synthesise an
/// alignment without an extension pass through here so the counts — and the
/// NM tag derived from them — describe the placement that was actually
/// written.
pub(crate) fn recount_mm_ea_in(
    reference: &crate::fasta::Reference,
    query: &[u8],
    aln: &mut ReadAlignment,
) {
    let contig = &reference.contigs[aln.contig as usize];
    let mut scratch: Vec<u8> = Vec::new();
    let (mut ro, mut rf) = (0u32, aln.pos);
    let (mut mm, mut ag, mut tc) = (0u32, 0u32, 0u32);
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) => {
                // Defensive: L1 projection can place a read past a
                // contig end (bad annotation); clamp instead of panic.
                let avail = (*n).min(contig.len.saturating_sub(rf));
                let s0 = scratch.len();
                if avail > 0 {
                    contig.decode_append(rf, rf + avail, &mut scratch);
                }
                let got = scratch.len() - s0;
                for k in 0..(*n as usize) {
                    let qb = query[ro as usize + k];
                    let rb = if k < got { scratch[s0 + k] } else { b'N' };
                    if qb != b'N' && rb != b'N' && qb != rb {
                        match crate::extend::exempt_dir(rb, qb) {
                            Some(true) => ag += 1,
                            Some(false) => tc += 1,
                            None => mm += 1,
                        }
                    }
                }
                ro += *n;
                rf += *n;
            }
            CigarOp::Ins(n) | CigarOp::SoftClip(n) => ro += *n,
            CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
        }
    }
    aln.mm_count = mm;
    aln.ea_ag = ag;
    aln.ea_tc = tc;
}

/** Reverse complement (N stays N).
 */
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    revcomp_into(seq, &mut out);
    out
}

/// [`revcomp`] into a reusable buffer (`out` is cleared, then filled).
pub fn revcomp_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(
        seq.iter()
            .rev()
            .map(|&b| crate::fasta::Base::from_ascii(b).complement().to_ascii()),
    );
}

/// Total read consumption (M/I/S) of a CIGAR.
fn cigar_read_span(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n,
            _ => 0,
        })
        .sum()
}

/// Tail-rescue: count base mismatches (+ indels) in the tail portion of a rescue CIGAR.
/// `right_tail = true` → tail after RefSkip; `false` → tail before RefSkip.
/// Non-zero = imperfect tail (any mismatch, insertion, or deletion).
fn count_tail_mismatches(
    cigar: &[CigarOp],
    pos: u32,
    contig: &crate::fasta::Contig,
    query: &[u8],
    right_tail: bool,
) -> u32 {
    let mut ref_pos = pos;
    let mut read_pos = 0u32;
    let mut seen_skip = false;
    let mut mismatches = 0u32;
    for op in cigar {
        let in_tail = if right_tail { seen_skip } else { !seen_skip };
        match op {
            CigarOp::RefSkip(n) => {
                seen_skip = true;
                ref_pos += n;
            }
            CigarOp::Match(n) => {
                if in_tail {
                    for i in 0..*n {
                        let qb = query[(read_pos + i) as usize];
                        // Out-of-bounds reference positions count as mismatches.
                        if ref_pos + i >= contig.len {
                            mismatches += 1;
                            continue;
                        }
                        let ra = contig.base(ref_pos + i).to_ascii();
                        if ra != b'N' && !qb.eq_ignore_ascii_case(&ra) {
                            mismatches += 1;
                        }
                    }
                }
                ref_pos += n;
                read_pos += n;
            }
            CigarOp::Ins(n) => {
                if in_tail {
                    mismatches += n;
                }
                read_pos += n;
            }
            CigarOp::Del(n) => {
                if in_tail {
                    mismatches += n;
                }
                ref_pos += n;
            }
            CigarOp::SoftClip(n) => {
                read_pos += n;
            }
        }
    }
    mismatches
}

/// Decoy contract: contigs whose name carries the `decoy_` prefix are
/// decoys. The prefix travels in the FASTA itself, so the index, the BAM
/// header and every downstream stage see it without a sidecar.
pub fn is_decoy_contig(reference: &crate::fasta::Reference, contig: u32) -> bool {
    reference
        .contigs
        .get(contig as usize)
        .is_some_and(|c| c.name.starts_with("decoy_"))
}

/// Reference-side M blocks (spliced gaps excluded) of a placement, as
/// sorted disjoint `(start, end)` intervals.
fn ref_m_blocks(pos: u32, cigar: &[CigarOp]) -> Vec<(u32, u32)> {
    let mut out: Vec<(u32, u32)> = Vec::new();
    let mut rf = pos;
    for op in cigar {
        match op {
            CigarOp::Match(n) => {
                out.push((rf, rf + *n));
                rf += *n;
            }
            CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
            _ => {}
        }
    }
    out
}

/// Footprint overlap of two placements: shared reference M bases divided by
/// the smaller footprint. Two candidates on the same contig and strand whose
/// overlap reaches 0.5 are the same locus; anything below is a genuine
/// alternative placement.
fn footprint_overlap(a: &Placement, b: &Placement) -> f64 {
    if a.contig != b.contig || a.strand != b.strand {
        return 0.0;
    }
    let ba = ref_m_blocks(a.pos, &a.cigar);
    let bb = ref_m_blocks(b.pos, &b.cigar);
    let la: u32 = ba.iter().map(|(s, e)| e - s).sum();
    let lb: u32 = bb.iter().map(|(s, e)| e - s).sum();
    if la == 0 || lb == 0 {
        return 0.0;
    }
    let mut inter: u32 = 0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < ba.len() && j < bb.len() {
        let lo = ba[i].0.max(bb[j].0);
        let hi = ba[i].1.min(bb[j].1);
        if hi > lo {
            inter += hi - lo;
        }
        if ba[i].1 < bb[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    f64::from(inter) / f64::from(la.min(lb))
}

/// Decision order: score desc, then fewer ordinary mismatches, more
/// informative query coverage, origin and coordinates (deterministic).
fn sort_placements(ps: &mut [Placement]) {
    ps.sort_by(|a, b| {
        b.strict
            .score
            .cmp(&a.strict.score)
            .then(a.strict.mm.cmp(&b.strict.mm))
            .then(
                b.informative_query_coverage()
                    .cmp(&a.informative_query_coverage()),
            )
            .then(a.origin.rank().cmp(&b.origin.rank()))
            .then(a.contig.cmp(&b.contig))
            .then(a.pos.cmp(&b.pos))
    });
}

/// The editing direction a placement's mismatch spectrum implies when it
/// looks hyperedited: single-direction edits at or above
/// [`crate::mapq::HYPEREDIT_MIN_DENSITY`] of the aligned columns with at most
/// [`crate::mapq::HYPEREDIT_MAX_NON_EDIT_MM`] ordinary mismatches (the same
/// criterion as `ReadAlignment::hyperedited_signature`). None otherwise.
fn hyper_direction(p: &Placement) -> Option<EditDir> {
    let s = &p.strict;
    if s.mm > crate::mapq::HYPEREDIT_MAX_NON_EDIT_MM || (s.ea_ag > 0 && s.ea_tc > 0) {
        return None;
    }
    let ea = s.ea_ag + s.ea_tc;
    if s.aligned == 0 || f64::from(ea) / f64::from(s.aligned) < crate::mapq::HYPEREDIT_MIN_DENSITY {
        return None;
    }
    Some(EditDir::from_plus(s.ea_ag > 0))
}

/// Keep the best-scored representative of each locus in each coverage class.
/// A short local explanation cannot erase a full-length fallback or its edit
/// signature; positive score is checked only after any hyper rescoring.
#[cfg(test)]
fn dedup_by_footprint(ps: &mut Vec<Placement>, read_len: usize) {
    dedup_by_footprint_with_options(ps, read_len, false);
}

fn dedup_by_footprint_with_options(
    ps: &mut Vec<Placement>,
    read_len: usize,
    exact_path_competition: bool,
) {
    sort_placements(ps);
    let mut kept: Vec<Placement> = Vec::with_capacity(ps.len());
    for p in ps.drain(..) {
        if kept.iter().all(|k| {
            k.coverage_eligible(read_len) != p.coverage_eligible(read_len)
                || footprint_overlap(k, &p) < 0.5
                || (exact_path_competition && exact_path::distinct(k, &p, read_len))
        }) {
            kept.push(p);
        }
    }
    *ps = kept;
}

mod exact_path;
mod short_body;

impl ReadAlignment {
    #[cfg(test)]
    fn mapq_of_checked(&self) -> u8 {
        crate::mapq::mapq_of(self)
    }
}

#[cfg(test)]
pub(crate) mod prefix_tests;

#[cfg(test)]
pub(crate) mod repeat_context_tests;

#[cfg(test)]
pub(crate) mod exact_path_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_from(seq: &str) -> crate::fasta::Reference {
        crate::fasta::parse_fasta_bytes(format!(">chrT\n{seq}\n").as_bytes()).expect("fasta")
    }

    fn aligner_on(seq: &str) -> Aligner<'static> {
        let reference = ref_from(seq);
        let index = Box::leak(Box::new(Index::build(
            reference,
            SeedParams { k: 15, w: 5 },
        )));
        Aligner::new(index, AlignConfig::rna_default())
    }

    fn placement(contig: u32, pos: u32, cigar: &[CigarOp]) -> Placement {
        Placement::from_genomic(
            contig,
            pos,
            Strand::Plus,
            cigar.to_vec(),
            Vec::new(),
            0,
            100,
        )
    }

    fn sparse_window_body(n: usize) -> Vec<u8> {
        let mut state = 0x75bcd15u64;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn sparse_window_chain(n: u32, near_head: bool, strand: Strand, pos: u32) -> Chain {
        let qs = if near_head {
            [2, 12, 27]
        } else {
            [n - 42, n - 32, n - 17]
        };
        Chain {
            contig: 0,
            strand,
            score: 40,
            anchors: qs
                .into_iter()
                .map(|qpos| crate::seed::Anchor {
                    contig: 0,
                    strand,
                    qpos,
                    rpos: pos + qpos,
                })
                .collect(),
        }
    }

    #[test]
    fn sparse_window_retains_full_ordinary_body() {
        for n in [150, 300] {
            for near_head in [true, false] {
                for strand in [Strand::Plus, Strand::Minus] {
                    let oriented = sparse_window_body(n);
                    let mut genome = vec![b'N'; 2000 + n];
                    genome[1000..1000 + n].copy_from_slice(&oriented);
                    let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                    let seq = if strand == Strand::Plus {
                        oriented.clone()
                    } else {
                        revcomp(&oriented)
                    };
                    let chain = sparse_window_chain(n as u32, near_head, strand, 1000);
                    let mut out = Vec::new();
                    let mut budget = EXTEND_BUDGET;
                    aligner.extend_unspliced_chain(
                        &seq,
                        &chain,
                        false,
                        &revcomp(&seq),
                        &mut budget,
                        &mut out,
                    );
                    assert_eq!(out.len(), 1, "{n} {near_head} {strand:?}");
                    aligner.rescore_strict(&oriented, &mut out[0]);
                    assert_eq!(
                        out[0].cigar,
                        vec![CigarOp::Match(n as u32)],
                        "{n} {near_head} {strand:?}"
                    );
                    assert_eq!(out[0].pos, 1000);
                    assert_eq!(out[0].strict.score, 2 * n as i32);
                    assert_eq!(out[0].strict.aligned, n as u32);
                    assert!(budget <= EXTEND_BUDGET);
                }
            }
        }
    }

    #[test]
    fn sparse_window_respects_contig_edges() {
        for n in [150, 300] {
            for (left, right) in [(0, 1000), (1000, 0)] {
                for near_head in [true, false] {
                    for strand in [Strand::Plus, Strand::Minus] {
                        let oriented = sparse_window_body(n);
                        let mut genome = vec![b'N'; left + n + right];
                        genome[left..left + n].copy_from_slice(&oriented);
                        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                        let seq = if strand == Strand::Plus {
                            oriented.clone()
                        } else {
                            revcomp(&oriented)
                        };
                        let chain = sparse_window_chain(n as u32, near_head, strand, left as u32);
                        let mut out = Vec::new();
                        let mut budget = EXTEND_BUDGET;
                        aligner.extend_unspliced_chain(
                            &seq,
                            &chain,
                            false,
                            &revcomp(&seq),
                            &mut budget,
                            &mut out,
                        );
                        assert_eq!(out.len(), 1);
                        assert_eq!(out[0].pos, left as u32);
                        assert_eq!(out[0].cigar, vec![CigarOp::Match(n as u32)]);
                    }
                }
            }
        }
    }

    #[test]
    fn sparse_window_negative_technical_tails_remain_clipped() {
        for n in [150, 300] {
            for tail_base in [b'N', b'C'] {
                for near_head in [true, false] {
                    for strand in [Strand::Plus, Strand::Minus] {
                        let body = sparse_window_body(n);
                        let oriented =
                            [vec![tail_base; 8], body.clone(), vec![tail_base; 8]].concat();
                        let m = oriented.len() as u32;
                        // Concrete A reference flanks make C tails ordinary
                        // negative-score mismatches, not ambiguous columns.
                        let mut genome = vec![b'A'; 2000 + n + 16];
                        genome[1008..1008 + n].copy_from_slice(&body);
                        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                        let seq = if strand == Strand::Plus {
                            oriented.clone()
                        } else {
                            revcomp(&oriented)
                        };
                        let mut chain = sparse_window_chain(m, near_head, strand, 1000);
                        let qs = if near_head {
                            [10, 20, 35]
                        } else {
                            [m - 50, m - 40, m - 25]
                        };
                        for (anchor, q) in chain.anchors.iter_mut().zip(qs) {
                            anchor.qpos = q;
                            anchor.rpos = 1000 + q;
                        }
                        let mut out = Vec::new();
                        let mut budget = EXTEND_BUDGET;
                        aligner.extend_unspliced_chain(
                            &seq,
                            &chain,
                            false,
                            &revcomp(&seq),
                            &mut budget,
                            &mut out,
                        );
                        assert_eq!(out.len(), 1);
                        assert_eq!(out[0].pos, 1008);
                        assert_eq!(
                            out[0].cigar,
                            vec![
                                CigarOp::SoftClip(8),
                                CigarOp::Match(n as u32),
                                CigarOp::SoftClip(8),
                            ]
                        );
                        aligner.rescore_strict(&oriented, &mut out[0]);
                        assert_eq!(out[0].strict.score, 2 * n as i32);
                    }
                }
            }
        }
    }

    #[test]
    fn sparse_window_preserves_neutral_terminal_edits_on_dp_path() {
        for n in [150, 300] {
            for dir in [EditDir::Ag, EditDir::Tc] {
                for near_head in [true, false] {
                    for strand in [Strand::Plus, Strand::Minus] {
                        let mut reference_body = sparse_window_body(n);
                        let (rb, qb) = if dir == EditDir::Ag {
                            (b'A', b'G')
                        } else {
                            (b'T', b'C')
                        };
                        reference_body[0] = rb;
                        reference_body[n - 1] = rb;
                        let mut oriented = reference_body.clone();
                        oriented[0] = qb;
                        oriented[n - 1] = qb;
                        // An internal ambiguous query base makes the short
                        // verification path use DP. Myers endpoint tie policy
                        // is a separate, pre-existing limitation.
                        if n <= 256 {
                            oriented[n / 2] = b'N';
                        }
                        let mut genome = vec![b'N'; 2000 + n];
                        genome[1000..1000 + n].copy_from_slice(&reference_body);
                        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                        aligner.config.extend.edit_dir = dir;
                        let seq = if strand == Strand::Plus {
                            oriented.clone()
                        } else {
                            revcomp(&oriented)
                        };
                        let chain = sparse_window_chain(n as u32, near_head, strand, 1000);
                        let mut out = Vec::new();
                        let mut budget = EXTEND_BUDGET;
                        aligner.extend_unspliced_chain(
                            &seq,
                            &chain,
                            false,
                            &revcomp(&seq),
                            &mut budget,
                            &mut out,
                        );
                        assert_eq!(out.len(), 1);
                        assert_eq!(out[0].pos, 1000);
                        assert_eq!(out[0].cigar, vec![CigarOp::Match(n as u32)]);
                        aligner.rescore_with(&oriented, &mut out[0], dir);
                        let uncounted = usize::from(n <= 256);
                        assert_eq!(out[0].strict.score, 2 * (n - 2 - uncounted) as i32);
                    }
                }
            }
        }
    }

    #[test]
    fn sparse_window_does_not_bypass_shared_dp_budget() {
        for n in [150, 300] {
            let seq = sparse_window_body(n);
            let mut genome = vec![b'N'; 2000 + n];
            genome[1000..1000 + n].copy_from_slice(&seq);
            let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
            let chain = sparse_window_chain(n as u32, true, Strand::Plus, 1000);
            let mut out = Vec::new();
            let mut budget = 0;
            aligner.extend_unspliced_chain(
                &seq,
                &chain,
                false,
                &revcomp(&seq),
                &mut budget,
                &mut out,
            );
            assert_eq!(budget, 0);
            if n <= 256 {
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].cigar, vec![CigarOp::Match(n as u32)]);
            } else {
                assert!(out.is_empty());
            }
        }
    }

    #[test]
    fn sparse_window_geometry_preserves_old_windows_and_bounds() {
        for n in [60u32, 150, 300, 1024] {
            let flank = n / 2 + 30;
            for qs in (0..=n - 15).step_by(7) {
                for qe in (qs + 15..=n).step_by(11) {
                    for pos in [0, 1, 1000] {
                        for extra in [0, 1000] {
                            let chain = Chain {
                                contig: 0,
                                strand: Strand::Plus,
                                score: 40,
                                anchors: [qs, qe - 15]
                                    .into_iter()
                                    .map(|qpos| crate::seed::Anchor {
                                        contig: 0,
                                        strand: Strand::Plus,
                                        qpos,
                                        rpos: pos + qpos,
                                    })
                                    .collect(),
                            };
                            let len = pos + n + extra;
                            let (start, end) = continuous_extension_window(&chain, 15, n, len);
                            assert!(start <= pos && end >= pos + n);
                            assert!(start <= end && end <= len);
                            assert!(end - start <= n + 2 * flank);
                            let (rs, re) = chain.ref_span(15);
                            if qs <= flank && n - qe <= flank {
                                assert_eq!(
                                    (start, end),
                                    (
                                        rs.saturating_sub(flank),
                                        (re + flank).min(len).min(rs + n + flank),
                                    )
                                );
                            }
                            // Displaced reference anchors cannot expand the
                            // continuous window to the entire reference gap.
                            let mut displaced = chain.clone();
                            displaced.anchors[1].rpos += 50000;
                            let (a, b) =
                                continuous_extension_window(&displaced, 15, n, len + 50000);
                            assert!(a <= b && b <= len + 50000);
                            assert!(b - a <= n + 2 * flank);
                        }
                    }
                }
            }
        }
        let mut chain = sparse_window_chain(150, true, Strand::Plus, 0);
        for anchor in &mut chain.anchors {
            anchor.rpos += u32::MAX - 100;
        }
        let (a, b) = continuous_extension_window(&chain, 15, 150, u32::MAX);
        assert!(a <= b && b - a <= 360);
    }

    // Synthetic rounds isolate scheduling from seed discovery. The complete
    // round has sparse anchors at the exact ordinary locus; a repeat round has
    // higher chain scores but poorer base evidence and includes the exact locus
    // after two competitors. Every short query contains N in its outer tail,
    // forcing verification to use the shared DP fallback. Long non-N queries
    // exercise the same budget without the <=256-base Myers path.
    fn round_competition(
        strand: Strand,
        leading: bool,
        long_non_n: bool,
    ) -> (Aligner<'static>, Vec<u8>, RoundOutcome, RoundOutcome, u32) {
        let mut state = 0x9e3779b9u64;
        let body_len = if long_non_n { 270 } else { 135 };
        let tail_len = if long_non_n { 30 } else { 15 };
        let body: Vec<u8> = (0..body_len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let tail = vec![if long_non_n { b'C' } else { b'N' }; tail_len];
        let oriented = if leading {
            [tail, body.clone()].concat()
        } else {
            [body.clone(), tail].concat()
        };
        let read = if strand == Strand::Plus {
            oriented
        } else {
            revcomp(&oriented)
        };
        let offset = if leading { tail_len as u32 } else { 0 };
        let strong_pos = 20000u32;
        let mut genome = vec![b'A'; 22000];
        genome[strong_pos as usize..strong_pos as usize + body_len].copy_from_slice(&body);
        let chain = |pos: u32, score: i32, qs: &[u32]| Chain {
            contig: 0,
            strand,
            score,
            anchors: qs
                .iter()
                .map(|&q| crate::seed::Anchor {
                    contig: 0,
                    strand,
                    qpos: offset + q,
                    rpos: pos + q,
                })
                .collect(),
        };
        let complete = RoundOutcome {
            chains: vec![chain(strong_pos, 60, &[45, 100])],
            second: 0,
        };
        let mut repeats = Vec::new();
        for i in 0..10 {
            let pos = 1000 + i * 1000;
            let target = &mut genome[pos as usize..pos as usize + body_len];
            target.copy_from_slice(&body);
            for q in [5, 55, 65, 115, 125] {
                let rb = target[q];
                target[q] = b"ACGT"
                    .iter()
                    .copied()
                    .find(|&b| b != rb && crate::extend::exempt_dir(b, rb).is_none())
                    .unwrap();
            }
            repeats.push(chain(pos, 120 - i as i32, &[20, 90]));
        }
        repeats.insert(4, chain(strong_pos, 116, &[20, 90]));
        let repeat = RoundOutcome {
            second: 119,
            chains: repeats,
        };
        (
            aligner_on(std::str::from_utf8(&genome).unwrap()),
            read,
            complete,
            repeat,
            strong_pos,
        )
    }

    #[test]
    fn limited_repeat_round_preserves_complete_path_before_shared_dp() {
        for strand in [Strand::Plus, Strand::Minus] {
            for leading in [false, true] {
                for long_non_n in [false, true] {
                    let (mut aligner, read, complete, repeat, strong_pos) =
                        round_competition(strand, leading, long_non_n);
                    let mins = crate::seed::minimizers(&read, aligner.config.seed);
                    let mut candidates = aligner.extend_genomic_rounds(
                        &read,
                        &mins,
                        Some(&complete),
                        Some(&repeat),
                        None,
                        true,
                    );
                    let query = if strand == Strand::Plus {
                        read.clone()
                    } else {
                        revcomp(&read)
                    };
                    for p in &mut candidates {
                        aligner.rescore_strict(&query, p);
                    }
                    assert!(
                        candidates.iter().any(|p| p.pos == strong_pos),
                        "complete path displaced before arbitration: {strand:?}, \
                         leading={leading}, long_non_n={long_non_n}"
                    );
                    assert_eq!(
                        candidates.len(),
                        EXTEND_BUDGET,
                        "all fixture paths require DP; the budget must remain two"
                    );
                    dedup_by_footprint(&mut candidates, read.len());
                    let budget = crate::index::AnchorBudget {
                        eligible_anchors: 25000,
                        retained_anchors: 19000,
                        omitted_seeds: 1,
                        repeat_context: None,
                    };
                    let result = aligner
                        .arbitrate_with_budget(&read, candidates, mins.len(), true, Some(budget))
                        .unwrap();
                    let (body, tail) = if long_non_n { (270, 30) } else { (135, 15) };
                    assert_eq!(result.pos, strong_pos);
                    assert_eq!(result.strand, strand);
                    assert_eq!(
                        result.score,
                        2 * body as i32,
                        "{strand:?} leading={leading} long_non_n={long_non_n} cigar={:?}",
                        result.cigar
                    );
                    let expected = if leading {
                        vec![CigarOp::SoftClip(tail), CigarOp::Match(body)]
                    } else {
                        vec![CigarOp::Match(body), CigarOp::SoftClip(tail)]
                    };
                    assert_eq!(
                        result.cigar, expected,
                        "outer technical tail must stay clipped"
                    );
                    assert_eq!(result.mapq_of_checked(), 0);
                }
            }
        }
    }

    #[test]
    fn complete_repeat_round_preserves_previous_extension_priority() {
        let (mut aligner, read, complete, repeat, strong_pos) =
            round_competition(Strand::Plus, false, false);
        let mins = crate::seed::minimizers(&read, aligner.config.seed);
        let candidates = aligner.extend_genomic_rounds(
            &read,
            &mins,
            Some(&complete),
            Some(&repeat),
            None,
            false,
        );
        assert_eq!(candidates.len(), EXTEND_BUDGET);
        assert!(
            candidates.iter().all(|p| p.pos != strong_pos),
            "scope control: complete search retains its existing scheduling"
        );
    }

    #[test]
    fn limited_repeat_round_can_beat_a_weaker_complete_path() {
        let (mut aligner, read, complete, repeat, strong_pos) =
            round_competition(Strand::Plus, true, false);
        let weak_complete = RoundOutcome {
            chains: vec![repeat.chains[0].clone()],
            second: 0,
        };
        let strong_repeat = RoundOutcome {
            chains: vec![complete.chains[0].clone()],
            second: 0,
        };
        let mins = crate::seed::minimizers(&read, aligner.config.seed);
        let mut candidates = aligner.extend_genomic_rounds(
            &read,
            &mins,
            Some(&weak_complete),
            Some(&strong_repeat),
            None,
            true,
        );
        assert!(candidates
            .iter()
            .any(|p| p.pos == weak_complete.chains[0].ref_span(15).0 - 20));
        assert!(candidates.iter().any(|p| p.pos == strong_pos));
        for p in &mut candidates {
            aligner.rescore_strict(&read, p);
        }
        dedup_by_footprint(&mut candidates, read.len());
        let result = aligner
            .arbitrate(&read, candidates, mins.len(), true)
            .unwrap();
        assert_eq!(
            result.pos, strong_pos,
            "strict evidence, not round origin, decides"
        );
        assert_eq!(result.score, 270);
    }

    #[test]
    fn limited_repeat_without_complete_chain_keeps_existing_search() {
        let (mut aligner, read, _, repeat, _) = round_competition(Strand::Minus, false, false);
        let mins = crate::seed::minimizers(&read, aligner.config.seed);
        let absent = aligner.extend_genomic_rounds(&read, &mins, None, Some(&repeat), None, true);
        let empty = RoundOutcome {
            chains: Vec::new(),
            second: 0,
        };
        let present_empty =
            aligner.extend_genomic_rounds(&read, &mins, Some(&empty), Some(&repeat), None, true);
        assert_eq!(absent.len(), EXTEND_BUDGET);
        assert_eq!(
            absent.iter().map(|p| (p.pos, &p.cigar)).collect::<Vec<_>>(),
            present_empty
                .iter()
                .map(|p| (p.pos, &p.cigar))
                .collect::<Vec<_>>()
        );
        assert!(aligner
            .extend_genomic_rounds(&read, &mins, None, None, None, true)
            .is_empty());
    }

    #[test]
    fn round_extension_keeps_eight_loci_and_exact_repeat_ambiguity() {
        let (fixture, _, _, _, _) = round_competition(Strand::Plus, false, false);
        let mut body = Vec::new();
        fixture.index.reference.contigs[0].decode_into(20000, 20135, &mut body);
        let mut genome = vec![b'N'; 14000];
        let chains: Vec<Chain> = (0..12)
            .map(|i| {
                let pos = 500 + i * 1000;
                genome[pos as usize..pos as usize + body.len()].copy_from_slice(&body);
                Chain {
                    contig: 0,
                    strand: Strand::Plus,
                    score: 135,
                    anchors: [20, 90]
                        .iter()
                        .map(|&q| crate::seed::Anchor {
                            contig: 0,
                            strand: Strand::Plus,
                            qpos: q,
                            rpos: pos + q,
                        })
                        .collect(),
                }
            })
            .collect();
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let complete = RoundOutcome {
            chains,
            second: 135,
        };
        let mins = crate::seed::minimizers(&body, aligner.config.seed);
        let mut candidates =
            aligner.extend_genomic_rounds(&body, &mins, Some(&complete), None, None, true);
        assert_eq!(candidates.len(), MAX_TIED_CANDIDATES);
        for p in &mut candidates {
            aligner.rescore_strict(&body, p);
        }
        dedup_by_footprint(&mut candidates, body.len());
        let result = aligner
            .arbitrate(&body, candidates, mins.len(), true)
            .unwrap();
        assert_eq!(result.n_cands, MAX_TIED_CANDIDATES);
        assert_eq!(
            result.second_score, result.score,
            "exact alternatives remain visible"
        );
        assert_eq!(result.mapq_of_checked(), 0);
    }

    #[test]
    fn limited_repeat_splice_can_beat_continuous_priority_round() {
        for strand in [Strand::Plus, Strand::Minus] {
            for leading in [false, true] {
                let mut state = 0x719ac57bu64;
                let mut genome: Vec<u8> = (0..6000)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        b"ACGT"[(state & 3) as usize]
                    })
                    .collect();
                genome[2064..2066].copy_from_slice(b"GT");
                genome[3998..4000].copy_from_slice(b"AG");
                genome[2063] = b'T';
                genome[4000] = b'A';
                let body = [&genome[2000..2064], &genome[4000..4071]].concat();
                genome[1000..1135].copy_from_slice(&body);
                for q in [40, 46, 52, 58, 64, 70, 76, 82] {
                    let qb = body[q];
                    genome[1000 + q] = b"ACGT"
                        .iter()
                        .copied()
                        .find(|&rb| rb != qb && crate::extend::exempt_dir(rb, qb).is_none())
                        .unwrap();
                }
                let offset = if leading { 15 } else { 0 };
                let query = if leading {
                    [b"N".repeat(15), body].concat()
                } else {
                    [body, b"N".repeat(15)].concat()
                };
                let read = if strand == Strand::Plus {
                    query.clone()
                } else {
                    revcomp(&query)
                };
                let make = |q, r| crate::seed::Anchor {
                    contig: 0,
                    strand,
                    qpos: offset + q,
                    rpos: r,
                };
                let complete = RoundOutcome {
                    chains: vec![Chain {
                        contig: 0,
                        strand,
                        score: 60,
                        anchors: vec![make(20, 1020), make(90, 1090)],
                    }],
                    second: 0,
                };
                let repeat = RoundOutcome {
                    chains: vec![Chain {
                        contig: 0,
                        strand,
                        score: 120,
                        anchors: [5, 25, 45, 70, 90, 110]
                            .iter()
                            .map(|&q| make(q, if q < 64 { 2000 + q } else { 4000 + q - 64 }))
                            .collect(),
                    }],
                    second: 0,
                };
                let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                let mins = crate::seed::minimizers(&read, aligner.config.seed);
                let mut candidates = aligner.extend_genomic_rounds(
                    &read,
                    &mins,
                    Some(&complete),
                    Some(&repeat),
                    None,
                    true,
                );
                assert!(
                    candidates.iter().any(|p| p.pos == 1000),
                    "priority proposal must remain"
                );
                for p in &mut candidates {
                    aligner.rescore_strict(&query, p);
                }
                dedup_by_footprint(&mut candidates, read.len());
                let result = aligner
                    .arbitrate(&read, candidates, mins.len(), true)
                    .unwrap();
                let mut expected = vec![
                    CigarOp::Match(64),
                    CigarOp::RefSkip(1936),
                    CigarOp::Match(71),
                ];
                if leading {
                    expected.insert(0, CigarOp::SoftClip(15));
                } else {
                    expected.push(CigarOp::SoftClip(15));
                }
                assert_eq!(
                    result.pos, 2000,
                    "stronger spliced second-round path must compete"
                );
                assert_eq!(result.cigar, expected);
            }
        }
    }

    #[test]
    fn repeat_budget_retains_repeat_paths_without_confidence_or_technical_tail() {
        let mut state = 0x71f3ab92u64;
        let body: Vec<u8> = (0..135)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let unit = [b"N".repeat(50), body.clone(), b"N".repeat(50)].concat();
        let genome = unit.repeat(700);
        let mut index = Index::build(
            ref_from(std::str::from_utf8(&genome).unwrap()),
            SeedParams { k: 15, w: 5 },
        );
        index.freq_cutoff = 1; // exact repeat seeds require the existing mid-occ round
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        for strand in [Strand::Plus, Strand::Minus] {
            for leading in [false, true] {
                let query = if leading {
                    [b"C".repeat(15), body.clone()].concat()
                } else {
                    [body.clone(), b"C".repeat(15)].concat()
                };
                let read = if strand == Strand::Plus {
                    query
                } else {
                    revcomp(&query)
                };
                let mins = crate::seed::minimizers(&read, aligner.config.seed);
                let (all, _) = crate::index::collect_anchors(&index, &mins, 150, 15, MID_OCC_CAP);
                assert!(all.len() > MID_ROUND_MAX_ANCHORS);
                let (result, cache) = aligner.align_read_full(&read);
                let result = result.expect("budgeted whole-seed search must retain repeat paths");
                let budget = cache.unwrap().seed_budget.unwrap();
                assert_eq!(budget.eligible_anchors, all.len());
                assert!(budget.retained_anchors <= MID_ROUND_MAX_ANCHORS);
                assert!(budget.limited());
                assert_eq!(result.strand, strand);
                assert_eq!(result.mapq_of_checked(), 0);
                let expected = if leading {
                    vec![CigarOp::SoftClip(15), CigarOp::Match(135)]
                } else {
                    vec![CigarOp::Match(135), CigarOp::SoftClip(15)]
                };
                assert_eq!(
                    result.cigar, expected,
                    "technical sequence must stay clipped"
                );
                assert_eq!(result.pos % unit.len() as u32, 50);
            }
        }
    }

    #[test]
    fn limited_seed_search_cannot_gain_confidence_from_a_single_path_or_mate() {
        let genome = "ACGT".repeat(500);
        let mut aligner = aligner_on(&genome);
        let read = genome.as_bytes()[100..250].to_vec();
        let mut p = placement(0, 100, &[CigarOp::Match(150)]);
        aligner.rescore_strict(&read, &mut p);
        let budget = crate::index::AnchorBudget {
            eligible_anchors: 25000,
            retained_anchors: 19000,
            omitted_seeds: 3,
            repeat_context: None,
        };
        let result = aligner
            .arbitrate_with_budget(&read, vec![p.clone()], 30, true, Some(budget))
            .unwrap();
        assert_eq!(result.mapq_of_checked(), 0);
        assert!(crate::evidence::evidence_tag(&result).contains("seed_search=limited"));
        let mut overridden = result;
        overridden.mapq_override = Some(60);
        assert_eq!(
            overridden.mapq_of_checked(),
            0,
            "later override cannot erase truncation"
        );

        // An empty independent pool must still carry search incompleteness
        // into the seed-free mate-window fallback and its evidence.
        let cache = CachedPlacements {
            placements: Vec::new(),
            n_seeds: 30,
            genomic_any: false,
            seed_budget: Some(budget),
        };
        let pool = cache.diagnostic_pool(aligner.index, 150, "independent");
        assert!(pool.seed_budget.unwrap().limited());
        let mut observed = None;
        let mate = aligner
            .align_read_with_mate_observed(&read, 0, 100, Some(&cache), Some(&mut observed))
            .unwrap();
        assert_eq!(mate.mapq_of_checked(), 0);
        assert!(mate.rescued);
        assert!(observed.unwrap().seed_budget.unwrap().limited());
    }

    #[test]
    fn diagnostic_96_98_99_preserves_evidence_without_promoting_either_mate() {
        let mut aligner = aligner_on(&"A".repeat(2000));
        for strand in [Strand::Plus, Strand::Minus] {
            for covered in [96, 98, 99] {
                let mut p = placement(
                    0,
                    100,
                    &[CigarOp::Match(covered), CigarOp::SoftClip(150 - covered)],
                );
                p.strand = strand;
                aligner.rescore_strict(&b"A".repeat(150), &mut p);
                let cache = CachedPlacements {
                    placements: vec![p],
                    n_seeds: 5,
                    genomic_any: true,
                    seed_budget: None,
                };
                let pool = cache.diagnostic_pool(aligner.index, 150, "independent");
                assert_eq!(pool.candidates.len(), 1);
                assert_eq!(pool.candidates[0].informative_query_bases, covered);
                assert_eq!(
                    pool.candidates[0].eligibility,
                    if covered == 99 {
                        "eligible"
                    } else {
                        "below_coverage"
                    }
                );
                let query = vec![if strand == Strand::Plus { b'A' } else { b'T' }; 150];
                let result = aligner.arbitrate(&query, cache.placements.clone(), 5, true);
                assert_eq!(result.is_some(), covered == 99);
            }
        }
    }

    #[test]
    fn paired_fragment_diagnostic_52_53_uses_effective_coverage_in_both_pools() {
        let mut aligner = aligner_on(&"A".repeat(2000));
        for strand in [Strand::Plus, Strand::Minus] {
            for covered in [52, 53] {
                let mut p = placement(
                    0,
                    100,
                    &[CigarOp::Match(covered), CigarOp::SoftClip(80 - covered)],
                );
                p.strand = strand;
                aligner.rescore_strict(&b"A".repeat(80), &mut p);
                let cache = CachedPlacements {
                    placements: vec![p],
                    n_seeds: 5,
                    genomic_any: true,
                    seed_budget: None,
                };
                let query = vec![if strand == Strand::Plus { b'A' } else { b'T' }; 80];
                let result = aligner.arbitrate(&query, cache.placements.clone(), 5, true);
                assert_eq!(result.is_some(), covered == 53);
                let pool = cache.diagnostic_pool(aligner.index, 80, "independent");
                let raw = crate::fastq::FastqRecord {
                    name: b"boundary".to_vec(),
                    seq: vec![b'A'; 150],
                    qual: vec![b'I'; 150],
                };
                let fragment = crate::fragment::FragmentEvidence {
                    fragment_length: 80,
                    overlap_bases: 80,
                    viable_offsets: 1,
                    adapter_family: "truseq",
                    adapter_bases_r1: 33,
                    adapter_bases_r2: 33,
                };
                let mut bytes = Vec::new();
                crate::diagnostics::write_row_with_fragment(
                    &mut bytes,
                    &raw,
                    Some(1),
                    result.is_some(),
                    pool.clone(),
                    Some(pool),
                    Some(fragment),
                )
                .unwrap();
                let row: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
                assert_eq!(row["schema_version"], 2);
                assert_eq!(row["read_length"], 150);
                assert_eq!(row["minimum_primary_query_bases"], 99);
                assert_eq!(row["effective_read_length"], 80);
                assert_eq!(row["effective_minimum_primary_query_bases"], 53);
                for name in ["independent", "mate_attempt"] {
                    assert_eq!(row[name]["scored_candidates"], 1);
                    assert_eq!(row[name]["excluded_nonpositive"], 0);
                    assert_eq!(row[name]["excluded_invalid_geometry"], 0);
                    let c = &row[name]["candidates"][0];
                    assert_eq!(c["informative_query_bases"], covered);
                    assert_eq!(c["common_score"], 2 * covered);
                    assert_eq!(
                        c["eligibility"],
                        if covered == 53 {
                            "eligible"
                        } else {
                            "below_coverage"
                        }
                    );
                    let expected = match (strand, covered) {
                        (Strand::Plus, 52) => "52M98S",
                        (Strand::Plus, _) => "53M97S",
                        (Strand::Minus, 52) => "70S52M28S",
                        (Strand::Minus, _) => "70S53M27S",
                    };
                    assert_eq!(c["cigar"], expected);
                }
            }
        }
    }

    #[test]
    fn diagnostic_geometry_scores_and_informative_coverage_are_not_primary_overrides() {
        let mut aligner = aligner_on(&"A".repeat(2000));
        let mut read = b"A".repeat(98);
        read.extend_from_slice(&b"N".repeat(52));
        let mut ambiguous = placement(0, 100, &[CigarOp::Match(150)]);
        aligner.rescore_strict(&read, &mut ambiguous);
        let mut deleted = placement(
            0,
            100,
            &[CigarOp::Match(98), CigarOp::Del(20), CigarOp::SoftClip(52)],
        );
        aligner.rescore_strict(&b"A".repeat(150), &mut deleted);
        let mut skipped = placement(
            0,
            100,
            &[
                CigarOp::Match(49),
                CigarOp::RefSkip(100),
                CigarOp::Match(49),
                CigarOp::SoftClip(52),
            ],
        );
        aligner.rescore_strict(&b"A".repeat(150), &mut skipped);
        let mut nonpositive = ambiguous.clone();
        nonpositive.strict.score = 0;
        let mut internal_clip = ambiguous.clone();
        internal_clip.cigar = vec![
            CigarOp::Match(70),
            CigarOp::SoftClip(10),
            CigarOp::Match(70),
        ];
        let mut out_of_bounds = ambiguous.clone();
        out_of_bounds.pos = 1950;
        let mut wrong_length = ambiguous.clone();
        wrong_length.cigar = vec![CigarOp::Match(149)];
        let mut insertion_only = ambiguous.clone();
        insertion_only.cigar = vec![CigarOp::Ins(150)];
        let cache = CachedPlacements {
            placements: vec![
                ambiguous,
                deleted,
                skipped,
                nonpositive,
                internal_clip,
                out_of_bounds,
                wrong_length,
                insertion_only,
            ],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let pool = cache.diagnostic_pool(aligner.index, 150, "independent");
        assert_eq!(pool.candidates.len(), 3);
        assert_eq!(pool.excluded_nonpositive, 1);
        assert_eq!(pool.excluded_invalid_geometry, 4);
        assert!(pool
            .candidates
            .iter()
            .all(|c| c.informative_query_bases == 98 && c.eligibility == "below_coverage"));
    }

    #[test]
    fn diagnostic_different_locus_short_leader_preserves_eligible_fallback_and_ambiguity() {
        let mut genome = b"A".repeat(2000);
        genome[1000..1020].fill(b'C');
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let query = b"A".repeat(150);
        let mut short = placement(0, 100, &[CigarOp::Match(98), CigarOp::SoftClip(52)]);
        let mut eligible = placement(0, 1000, &[CigarOp::Match(150)]);
        aligner.rescore_strict(&query, &mut short);
        aligner.rescore_strict(&query, &mut eligible);
        let mut placements = vec![eligible, short];
        sort_placements(&mut placements);
        let cache = CachedPlacements {
            placements: placements.clone(),
            n_seeds: 10,
            genomic_any: true,
            seed_budget: None,
        };
        let before = aligner
            .arbitrate(&query, placements.clone(), 10, true)
            .unwrap();
        let pool = cache.diagnostic_pool(aligner.index, 150, "independent");
        assert_eq!(pool.candidates[0].eligibility, "below_coverage");
        assert_eq!(pool.candidates[1].eligibility, "eligible");
        let after = aligner.arbitrate(&query, placements, 10, true).unwrap();
        assert_eq!(
            (
                before.pos,
                before.score,
                before.second_score,
                before.mapq_of_checked()
            ),
            (
                after.pos,
                after.score,
                after.second_score,
                after.mapq_of_checked()
            )
        );
        assert_eq!(after.pos, 1000);
        assert_eq!(after.second_score, 196);
        assert_eq!(after.mapq_of_checked(), 0);
    }

    #[test]
    fn right_split_rescue_updates_position_and_short_tail_gate() {
        let j = Junction {
            contig: 0,
            start: 60,
            end: 160,
            minus_strand: false,
        };
        for tail_len in [10, 40] {
            let seq = format!(
                "{}{}GT{}AG{}{}",
                "A".repeat(20),
                "C".repeat(40),
                "T".repeat(96),
                "G".repeat(tail_len),
                "T".repeat(60 - tail_len)
            );
            let mut aligner = aligner_on(&seq).with_lib(Arc::new(JunctionLib::build([j, j])));
            let read = format!("{}{}", "C".repeat(40), "G".repeat(tail_len));
            // Reextension shifts the body start from 23 to 20. The short-tail
            // mismatch gate must also use that new start.
            let aln = ReadAlignment {
                pos: 23,
                cigar: vec![CigarOp::Match(40), CigarOp::SoftClip(tail_len as u32)],
                ..Default::default()
            };
            let rescued = aligner.apply_split_rescue(read.as_bytes(), aln);
            assert_eq!(rescued.pos, 20, "tail length {tail_len}: {rescued:?}");
            assert_eq!(
                rescued.cigar,
                vec![
                    CigarOp::Match(40),
                    CigarOp::RefSkip(100),
                    CigarOp::Match(tail_len as u32)
                ]
            );
            assert_eq!(rescued.junctions.len(), 1);
        }
    }

    #[test]
    fn rejected_split_keeps_original_primary() {
        let seq = format!(
            "{}{}GT{}AG{}",
            "A".repeat(20),
            "C".repeat(40),
            "T".repeat(96),
            "G".repeat(60)
        );
        let j = Junction {
            contig: 0,
            start: 60,
            end: 160,
            minus_strand: false,
        };
        let mut aligner = aligner_on(&seq).with_lib(Arc::new(JunctionLib::build([j, j])));
        let read = format!("{}AAAA{}", "C".repeat(40), "G".repeat(40));
        let original = vec![CigarOp::Match(40), CigarOp::SoftClip(44)];
        let aln = ReadAlignment {
            pos: 20,
            cigar: original.clone(),
            score: 80,
            ..Default::default()
        };
        let result = aligner.apply_split_rescue(read.as_bytes(), aln);
        assert_eq!(result.pos, 20);
        assert_eq!(result.cigar, original);
        assert_eq!(result.score, 80);
        assert!(result.junctions.is_empty());
    }

    #[test]
    fn marginal_short_tail_rescue_does_not_replace_the_primary() {
        // AT..AC intron of 22kb. Moving two already matched body bases
        // permits six clipped bases to match, but the existing signal/length
        // score improves by only three: below the split rewrite cost.
        let seq = format!(
            "{}{}AT{}ACATGGGGGG{}",
            "A".repeat(20),
            "C".repeat(38),
            "T".repeat(21996),
            "A".repeat(60)
        );
        let read = format!("{}ATGGGGGGCC", "C".repeat(38));
        let mut aligner = aligner_on(&seq);
        let original = vec![CigarOp::Match(40), CigarOp::SoftClip(8)];
        let aln = ReadAlignment {
            pos: 20,
            score: 80,
            cigar: original.clone(),
            ..Default::default()
        };
        let result = aligner.apply_split_rescue(read.as_bytes(), aln);
        assert_eq!(
            result.cigar, original,
            "weak tail must not create a discovery: {result:?}"
        );
        assert!(result.junctions.is_empty());
    }

    #[test]
    fn span_rescue_requires_sequence_gain_beyond_the_junction_prior() {
        // Both the contiguous GT and the acceptor GT exactly explain the
        // same two bases. Library support is a prior, not extra read evidence.
        let seq = format!(
            "{}GT{}AGGT{}",
            "C".repeat(60),
            "T".repeat(96),
            "A".repeat(60)
        );
        let j = Junction {
            contig: 0,
            start: 60,
            end: 160,
            minus_strand: false,
        };
        let mut aligner = aligner_on(&seq).with_lib(Arc::new(JunctionLib::build([j, j])));
        let read = format!("{}GT", "C".repeat(40));
        let original = vec![CigarOp::Match(42)];
        let aln = ReadAlignment {
            pos: 20,
            score: 84,
            cigar: original.clone(),
            ..Default::default()
        };
        let result = aligner.apply_span_rescue(read.as_bytes(), aln);
        assert_eq!(result.cigar, original);
        assert!(result.junctions.is_empty());
    }

    #[test]
    fn supported_span_rescue_keeps_a_genuine_sequence_improvement() {
        let seq = format!(
            "{}GT{}AGGG{}",
            "C".repeat(60),
            "T".repeat(96),
            "A".repeat(60)
        );
        let j = Junction {
            contig: 0,
            start: 60,
            end: 160,
            minus_strand: false,
        };
        let mut aligner = aligner_on(&seq).with_lib(Arc::new(JunctionLib::build([j, j])));
        let read = format!("{}GG", "C".repeat(40));
        let aln = ReadAlignment {
            pos: 20,
            cigar: vec![CigarOp::Match(42)],
            ..Default::default()
        };
        let result = aligner.apply_span_rescue(read.as_bytes(), aln);
        assert_eq!(
            result.cigar,
            vec![CigarOp::Match(40), CigarOp::RefSkip(100), CigarOp::Match(2)]
        );
        assert_eq!(result.junctions.len(), 1);
    }

    #[test]
    fn rescue_margin_uses_strict_scores_in_both_read_orientations() {
        let seq = format!("{}GC{}AG{}", "C".repeat(60), "T".repeat(96), "G".repeat(60));
        for strand in [Strand::Plus, Strand::Minus] {
            for tail in [5, 6] {
                let query = format!("{}{}", "C".repeat(40), "G".repeat(tail)).into_bytes();
                let read = if strand == Strand::Minus {
                    revcomp(&query)
                } else {
                    query
                };
                let mut aligner = aligner_on(&seq);
                let primary = ReadAlignment {
                    pos: 20,
                    strand,
                    score: 80,
                    cigar: vec![CigarOp::Match(40), CigarOp::SoftClip(tail as u32)],
                    ..Default::default()
                };
                let candidate = ReadAlignment {
                    pos: 20,
                    strand,
                    score: 10000, // arbitrary producer domain
                    cigar: vec![
                        CigarOp::Match(40),
                        CigarOp::RefSkip(100),
                        CigarOp::Match(tail as u32),
                    ],
                    junctions: vec![RefinedJunction {
                        junction: Junction {
                            contig: 0,
                            start: 60,
                            end: 160,
                            minus_strand: strand == Strand::Minus,
                        },
                        signal: SpliceSignal::GcAg,
                        known_support: 0,
                    }],
                    ..Default::default()
                };
                let result = aligner.select_rescue(&read, primary, candidate);
                if tail == 5 {
                    // 90 base score + 2 motif - 8 length = 84. The four-point
                    // improvement only ties the existing rewrite cost.
                    assert!(result.junctions.is_empty());
                    assert_eq!(result.score, 80);
                } else {
                    assert_eq!(result.junctions.len(), 1);
                    assert_eq!(result.score, 86);
                }
            }
        }
    }

    #[test]
    fn segmented_chain_keeps_the_contiguous_interpretation_of_each_exon() {
        // A repeat-derived first anchor joins a much stronger downstream
        // segment. Splice construction must not hide a full-length contiguous
        // placement at that segment; strict arbitration, not generation
        // order or mate fallback, decides which explanation wins.
        let mut state = 0x6a09e667u64;
        let mut random = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                })
                .collect()
        };
        let mut genome = random(16000);
        let mut read = random(150);
        read[47..51].copy_from_slice(b"AGGC");
        genome[12000..12150].copy_from_slice(&read);
        genome[2000..2049].copy_from_slice(&read[..49]);
        genome[2000..2008].fill(b'N');
        genome[2012] = b"ACGT".iter().copied().find(|&b| b != read[12]).unwrap();
        genome[2023] = b"ACGT".iter().copied().find(|&b| b != read[23]).unwrap();
        genome[2049..2051].copy_from_slice(b"GC");
        // One A>G substitution at the correct locus.
        genome[12030] = b'A';
        read[30] = b'G';
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut anchors = vec![crate::seed::Anchor {
            qpos: 35,
            rpos: 2035,
            contig: 0,
            strand: Strand::Plus,
        }];
        for q in [53, 60, 83, 95, 106, 118] {
            anchors.push(crate::seed::Anchor {
                qpos: q,
                rpos: 12000 + q,
                contig: 0,
                strand: Strand::Plus,
            });
        }
        let chain = Chain {
            contig: 0,
            strand: Strand::Plus,
            anchors,
            score: 90,
        };
        let mut candidates = Vec::new();
        aligner.extend_chain(
            &read,
            &crate::seed::minimizers(&read, aligner.config.seed),
            &chain,
            true,
            false,
            &revcomp(&read),
            &mut 8,
            false,
            &mut candidates,
        );
        assert!(
            candidates
                .iter()
                .any(|p| p.cigar.iter().any(|op| matches!(op, CigarOp::RefSkip(_)))),
            "fixture must also produce a spliced interpretation"
        );
        assert!(
            candidates
                .iter()
                .any(|p| p.pos == 12000 && p.cigar == vec![CigarOp::Match(150)]),
            "successful splicing must not suppress the contiguous candidate"
        );
        for p in &mut candidates {
            aligner.rescore_strict(&read, p);
        }
        dedup_by_footprint(&mut candidates, read.len());
        sort_placements(&mut candidates);
        assert_eq!(candidates[0].pos, 12000);
        assert_eq!(candidates[0].cigar, vec![CigarOp::Match(150)]);
        assert_eq!(candidates[0].strict.score, 294);
    }

    #[test]
    fn reused_buffers_two_reads_independent() {
        // One Aligner (one set of reusable buffers), consecutive reads with
        // different content: no state may leak between calls, and a replay
        // of the first read must reproduce the first result exactly
        // (allocation-elimination invariant).
        let contig = "GATTACAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGCAAGTGGCATTGACCAAGTTTGCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATT";
        let mut aln = aligner_on(contig);
        let r1 = &contig.as_bytes()[0..60];
        let r2 = &contig.as_bytes()[91..151];

        let a1 = aln.align_read(r1).expect("first read placed");
        assert_eq!(a1.pos, 0);

        let mut r2n = r2.to_vec();
        r2n[30] = b'N'; // different content through the same scratch buffers
        let a2 = aln.align_read(&r2n).expect("second read placed");
        assert_eq!(a2.pos, 91);
        assert_ne!(a1.pos, a2.pos);

        // Replay of the first read: identical output to the first call.
        let a1b = aln.align_read(r1).expect("replay placed");
        let mq1b = a1b.mapq_of_checked();
        let mq1 = a1.mapq_of_checked();
        assert_eq!(
            (a1b.pos, a1b.strand, a1b.score, a1b.cigar),
            (a1.pos, a1.strand, a1.score, a1.cigar)
        );
        assert_eq!(mq1b, mq1);
    }

    #[test]
    fn rescore_strict_manual_cigars() {
        let mut aln = aligner_on(&"ACGT".repeat(50));
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        // 50M perfect at pos 0: 50 matches.
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(read, &mut p);
        assert_eq!(p.strict.score, 100);
        assert_eq!(p.strict.matches, 50);
        assert_eq!(p.strict.mm, 0);
        assert_eq!(p.strict.aligned, 50);
        // One ordinary mismatch: 49*2 - 4.
        let mut r2 = read.to_vec();
        r2[25] = if r2[25] == b'A' { b'C' } else { b'A' };
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(&r2, &mut p);
        assert_eq!(p.strict.score, 49 * 2 - 4);
        assert_eq!(p.strict.mm, 1);
        // An A>G editing-type substitution is charged identically but
        // counted in the editing channel.
        let mut r3 = read.to_vec();
        let ref_b = aln.index.reference.contigs[0].base(25).to_ascii();
        if ref_b == b'A' {
            r3[25] = b'G'; // ref A / read G
        } else {
            // shift to a position with ref A: "ACGT"*n has A at multiples of 4
            let pos = 24; // A at 24 (0-based) in "ACGT" repeats
            let _ = pos;
            r3[24] = b'G';
        }
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(&r3, &mut p);
        assert_eq!(p.strict.score, 49 * 2 - 4);
        assert_eq!(p.strict.mm, 0);
        assert!(p.strict.ea_ag + p.strict.ea_tc == 1);
        // Insertion of 2: 48 aligned matches, -(4 + 2*2).
        let mut ins_read: Vec<u8> = read[..24].to_vec();
        ins_read.extend_from_slice(b"TT");
        ins_read.extend_from_slice(&read[24..48]);
        let mut p = placement(
            0,
            0,
            &[CigarOp::Match(24), CigarOp::Ins(2), CigarOp::Match(24)],
        );
        aln.rescore_strict(&ins_read, &mut p);
        assert_eq!(p.strict.score, 48 * 2 - (4 + 2 * 2));
        assert_eq!(p.strict.ins, 1);
    }

    #[test]
    fn rescore_strict_splice_signal() {
        // Exon1 | intron(GT + 56N + AG) | exon2; read spans 20+20 exon bases.
        let exon1 = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCAAGTTGACGTTACAGGCATTAC"; // 57
        let intron = format!("GT{}AG", "N".repeat(56)); // 60 total, GT..AG inside
        let exon2 = "CATGGACTTGAACCTTGACGTACATGGACTTGAACCTTGACGTA"; // 44
        let seq = format!("{exon1}{intron}{exon2}");
        let mut aln = aligner_on(&seq);
        let start = (exon1.len() - 20) as u32; // last 20 exon1 bases
        let read: Vec<u8> = {
            let s = seq.bytes().collect::<Vec<u8>>();
            let mut r: Vec<u8> = s[(start as usize)..exon1.len()].to_vec();
            let e2 = exon1.len() + intron.len();
            r.extend_from_slice(&s[e2..e2 + 20]);
            r
        };
        let n = intron.len() as u32;
        let mut p = placement(
            0,
            start,
            &[CigarOp::Match(20), CigarOp::RefSkip(n), CigarOp::Match(20)],
        );
        aln.rescore_strict(&read, &mut p);
        let expected = 40 * 2 + SpliceParams::default().gt_ag_bonus - intron_length_penalty(n);
        assert_eq!(p.strict.score, expected);
        assert_eq!(p.strict.n_splice, 1);
        assert_eq!(p.strict.aligned, 40);
        // A read antisense to the transcript crosses the same GT..AG intron
        // (seen as CT..AC from its side); the strict score must not depend
        // on which strand the read was sequenced from.
        let read_rc = revcomp(&read);
        let mut pm = placement(
            0,
            start,
            &[CigarOp::Match(20), CigarOp::RefSkip(n), CigarOp::Match(20)],
        );
        pm.strand = Strand::Minus;
        aln.rescore_strict(&revcomp(&read_rc), &mut pm);
        assert_eq!(pm.strict.score, expected);
    }

    #[test]
    fn exact_kmer_support_counts() {
        let mut aln = aligner_on(&"ACGT".repeat(50));
        let read: Vec<u8> = b"ACGT".repeat(13)[..50].to_vec();
        let p = placement(0, 0, &[CigarOp::Match(50)]);
        assert_eq!(aln.exact_kmer_support(&read, &p), 3); // 0,15,30
                                                          // One mismatch in the second 15-mer window pushes the scan past it.
        let mut r2 = read.clone();
        r2[20] = if r2[20] == b'A' { b'C' } else { b'A' };
        let n = aln.exact_kmer_support(&r2, &p);
        assert!(n < 3);
    }

    #[test]
    fn footprint_overlap_boundaries() {
        let a = placement(0, 0, &[CigarOp::Match(100)]);
        let half = placement(0, 50, &[CigarOp::Match(100)]);
        assert!(footprint_overlap(&a, &half) >= 0.5); // exactly 0.5: same locus
        let below = placement(0, 51, &[CigarOp::Match(100)]);
        assert!(footprint_overlap(&a, &below) < 0.5);
        let disjoint = placement(0, 200, &[CigarOp::Match(100)]);
        assert_eq!(footprint_overlap(&a, &disjoint), 0.0);
        // Different strand never overlaps.
        let mut other_strand = placement(0, 0, &[CigarOp::Match(100)]);
        other_strand.strand = Strand::Minus;
        assert_eq!(footprint_overlap(&a, &other_strand), 0.0);
        // Spliced placements overlap on M blocks only.
        let spliced = placement(
            0,
            0,
            &[
                CigarOp::Match(50),
                CigarOp::RefSkip(1000),
                CigarOp::Match(50),
            ],
        );
        let head = placement(0, 0, &[CigarOp::Match(50)]);
        assert!(footprint_overlap(&spliced, &head) > 0.9);
    }

    #[test]
    fn dedup_keeps_best_strict_per_locus() {
        let mut a = placement(0, 0, &[CigarOp::Match(100)]);
        a.strict.score = 150;
        let mut b = placement(0, 10, &[CigarOp::Match(100)]);
        b.strict.score = 200; // same locus, better
        let mut c = placement(0, 5000, &[CigarOp::Match(100)]);
        c.strict.score = 100; // different locus
        let mut ps = vec![a, b.clone(), c];
        dedup_by_footprint(&mut ps, 100);
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].strict.score, 200);
        assert_eq!(ps[0].pos, b.pos);
        assert_eq!(ps[1].strict.score, 100);
    }

    #[test]
    fn coverage_short_same_locus_cannot_erase_full_fallback() {
        let mut genome = b"A".repeat(1000);
        genome[..30].fill(b'C');
        let read = b"A".repeat(150);
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut full = placement(0, 0, &[CigarOp::Match(150)]);
        let mut short = placement(0, 70, &[CigarOp::SoftClip(70), CigarOp::Match(80)]);
        aligner.rescore_strict(&read, &mut full);
        aligner.rescore_strict(&read, &mut short);
        assert_eq!((full.strict.score, short.strict.score), (120, 160));
        assert_eq!(footprint_overlap(&full, &short), 1.0);
        let mut candidates = vec![short, full];
        dedup_by_footprint(&mut candidates, read.len());
        let result = aligner
            .arbitrate(&read, candidates, 10, true)
            .expect("an ineligible short hit must not delete the full placement");
        assert_eq!(result.cigar, vec![CigarOp::Match(150)]);
        assert_eq!(result.score, 120);
        assert_eq!(result.n_cands, 1);
        assert!(result.mapq_of_checked() < crate::mapq::PILEUP_GATE);
    }

    #[test]
    fn coverage_disjoint_short_hit_allows_fallback_but_retains_ambiguity() {
        let mut genome = b"A".repeat(1000);
        genome[500..530].fill(b'C');
        let read = b"A".repeat(150);
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut full = placement(0, 500, &[CigarOp::Match(150)]);
        let mut short = placement(0, 0, &[CigarOp::Match(80), CigarOp::SoftClip(70)]);
        aligner.rescore_strict(&read, &mut full);
        aligner.rescore_strict(&read, &mut short);
        let mut candidates = vec![full, short];
        dedup_by_footprint(&mut candidates, read.len());
        let result = aligner
            .arbitrate(&read, candidates, 10, true)
            .expect("continue past a disjoint coverage-ineligible score leader");
        assert_eq!(result.pos, 500);
        assert_eq!(result.score, 120);
        assert_eq!(result.second_score, 160);
        assert_eq!(result.mapq_of_checked(), 0);
        assert_eq!(result.n_cands, 1);
    }

    #[test]
    fn coverage_full_nonpositive_candidate_is_not_a_fallback() {
        let mut aligner = aligner_on(&"C".repeat(1000));
        let read = b"A".repeat(150);
        let mut full = placement(0, 0, &[CigarOp::Match(150)]);
        aligner.rescore_strict(&read, &mut full);
        assert_eq!(full.strict.score, -600);
        assert!(aligner
            .arbitrate(&read, vec![full.clone()], 10, true)
            .is_none());
        full.strict.score = 0;
        assert!(aligner.arbitrate(&read, vec![full], 10, true).is_none());
    }

    #[test]
    fn coverage_boundary_counts_query_evidence_in_both_orientations() {
        let mut aligner = aligner_on(&"A".repeat(1000));
        for strand in [Strand::Plus, Strand::Minus] {
            for covered in [98, 99] {
                let mut p = placement(
                    0,
                    0,
                    &[CigarOp::Match(covered), CigarOp::SoftClip(150 - covered)],
                );
                p.strand = strand;
                aligner.rescore_strict(&b"A".repeat(150), &mut p);
                let read = if strand == Strand::Plus { b'A' } else { b'T' };
                let result = aligner.arbitrate(&vec![read; 150], vec![p], 10, true);
                assert_eq!(result.is_some(), covered == 99);
            }
        }
        let mut with_insertion = placement(
            0,
            0,
            &[CigarOp::Match(90), CigarOp::Ins(9), CigarOp::SoftClip(51)],
        );
        aligner.rescore_strict(&b"A".repeat(150), &mut with_insertion);
        assert!(aligner
            .arbitrate(&b"A".repeat(150), vec![with_insertion], 10, true)
            .is_some());
    }

    #[test]
    fn coverage_deletions_skips_and_ambiguous_columns_do_not_fill_the_read() {
        let mut aligner = aligner_on(&"A".repeat(1000));
        let read = b"A".repeat(150);
        for gap in [CigarOp::Del(50), CigarOp::RefSkip(50)] {
            let mut p = placement(0, 0, &[CigarOp::Match(98), gap, CigarOp::SoftClip(52)]);
            aligner.rescore_strict(&read, &mut p);
            assert!(p.strict.score > 0);
            assert!(
                aligner.arbitrate(&read, vec![p], 10, true).is_none(),
                "reference gaps cannot satisfy query coverage"
            );
        }
        let mut ambiguous_read = b"A".repeat(98);
        ambiguous_read.extend_from_slice(&b"N".repeat(52));
        let mut p = placement(0, 0, &[CigarOp::Match(150)]);
        aligner.rescore_strict(&ambiguous_read, &mut p);
        assert!(p.strict.score > 0);
        assert!(aligner
            .arbitrate(&ambiguous_read, vec![p], 10, true)
            .is_none());
    }

    #[test]
    fn coverage_hyper_rescoring_can_promote_a_retained_full_candidate() {
        let mut aligner = aligner_on(&"A".repeat(1000));
        let mut read = b"A".repeat(80);
        read.extend_from_slice(&b"G".repeat(70));
        let mut full = placement(0, 0, &[CigarOp::Match(150)]);
        let mut short = placement(0, 0, &[CigarOp::Match(80), CigarOp::SoftClip(70)]);
        aligner.rescore_strict(&read, &mut full);
        aligner.rescore_strict(&read, &mut short);
        assert!(full.strict.score < 0 && hyper_direction(&full).is_some());
        let mut candidates = vec![full, short];
        dedup_by_footprint(&mut candidates, read.len());
        aligner.hyper_rearbitrate(&read, &revcomp(&read), &mut candidates);
        let result = aligner
            .arbitrate(&read, candidates, 10, true)
            .expect("retain a full candidate as hyperediting evidence before dedup");
        assert_eq!(result.cigar, vec![CigarOp::Match(150)]);
        assert_eq!(result.ea_ag, 70);
    }

    #[test]
    fn coverage_mate_window_without_eligible_candidate_returns_none() {
        let mut aligner = aligner_on(&"A".repeat(1000));
        let mut read = b"A".repeat(98);
        read.extend_from_slice(&b"G".repeat(52));
        let mut win = aligner
            .mate_window_placement(&read, 0, 100)
            .expect("reduced-alphabet window search supplies a partial candidate");
        aligner.rescore_strict(&read, &mut win);
        assert!(!win.primary_eligible(read.len()));
        let cache = CachedPlacements {
            placements: Vec::new(),
            n_seeds: 0,
            genomic_any: false,
            seed_budget: None,
        };
        assert!(aligner
            .align_read_with_mate(&read, 0, 100, Some(&cache))
            .is_none());
        let mut evidence = None;
        assert!(aligner
            .align_read_with_mate_observed(&read, 0, 100, Some(&cache), Some(&mut evidence))
            .is_none());
        let pool = evidence.unwrap();
        assert_eq!(pool.stage, "mate_window_rejected");
        assert!(!pool.candidates.is_empty());
        assert!(pool
            .candidates
            .iter()
            .all(|p| p.eligibility == "below_coverage"));
        let good = aligner
            .align_read_with_mate(&b"A".repeat(150), 0, 100, Some(&cache))
            .expect("eligible mate-window placements still work");
        assert!(good.rescued);
        assert_eq!(good.cigar, vec![CigarOp::Match(150)]);
        assert!(good.mapq_of_checked() <= crate::mapq::PILEUP_GATE);
        let mut evidence = None;
        let observed = aligner
            .align_read_with_mate_observed(
                &b"A".repeat(150),
                0,
                100,
                Some(&cache),
                Some(&mut evidence),
            )
            .unwrap();
        assert_eq!(observed.cigar, good.cigar);
        assert_eq!(observed.pos, good.pos);
        assert_eq!(observed.mapq_of_checked(), good.mapq_of_checked());
        assert_eq!(evidence.unwrap().stage, "mate_window_arbitration");
        let mut evidence = None;
        assert!(aligner
            .align_read_with_mate_observed(&[], 0, 100, Some(&cache), Some(&mut evidence))
            .is_none());
        let pool = evidence.unwrap();
        assert_eq!(pool.stage, "mate_window_no_proposal");
        assert!(pool.candidates.is_empty());
    }

    #[test]
    fn coverage_dedup_fallback_is_independent_of_generator_order() {
        let mut genome = b"A".repeat(1000);
        genome[..30].fill(b'C');
        let read = b"A".repeat(150);
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut full = placement(0, 0, &[CigarOp::Match(150)]);
        let mut short = placement(0, 70, &[CigarOp::SoftClip(70), CigarOp::Match(80)]);
        let mut alternative = placement(0, 500, &[CigarOp::Match(99), CigarOp::SoftClip(51)]);
        aligner.rescore_strict(&read, &mut full);
        aligner.rescore_strict(&read, &mut short);
        aligner.rescore_strict(&read, &mut alternative);
        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let pool = [full.clone(), short.clone(), alternative.clone()];
            let mut candidates = order.map(|i| pool[i].clone()).to_vec();
            dedup_by_footprint(&mut candidates, read.len());
            let result = aligner
                .arbitrate(&read, candidates, 10, true)
                .expect("eligible");
            assert_eq!(result.pos, 500);
            assert_eq!(result.second_score, 160);
            assert_eq!(result.n_cands, 2);
        }
    }

    #[test]
    fn neutral_same_locus_ties_keep_more_query_evidence_in_any_order() {
        let mut genome = b"A".repeat(1000);
        genome[649] = b'T';
        let mut query = b"A".repeat(150);
        query[149] = b'C';
        for strand in [Strand::Plus, Strand::Minus] {
            let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
            let read = if strand == Strand::Plus {
                query.clone()
            } else {
                revcomp(&query)
            };
            let mut full = placement(0, 500, &[CigarOp::Match(150)]);
            let mut clipped = placement(0, 500, &[CigarOp::Match(147), CigarOp::SoftClip(3)]);
            full.strand = strand;
            clipped.strand = strand;
            aligner.rescore_strict(&query, &mut full);
            aligner.rescore_strict(&query, &mut clipped);
            assert_eq!(full.strict.score, clipped.strict.score);
            assert_eq!(full.strict.mm, clipped.strict.mm);
            for mut candidates in [
                vec![clipped.clone(), full.clone()],
                vec![full.clone(), clipped.clone()],
            ] {
                dedup_by_footprint(&mut candidates, read.len());
                let result = aligner.arbitrate(&read, candidates, 10, true).unwrap();
                assert_eq!(
                    result.cigar, full.cigar,
                    "a same-score clipped tie must not erase terminal edit evidence"
                );
                assert_eq!(result.ea_tc, 1);
            }
        }
    }

    #[test]
    fn splice_competition_retains_neutral_terminal_edit_on_both_mates() {
        let mut state = 0x517cc1b7u64;
        let mut genome: Vec<u8> = (0..4000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        genome[1078..1080].copy_from_slice(b"CT");
        genome[2998..3000].copy_from_slice(b"AC");
        genome[1002] = b'T';
        let mut query = b"N".repeat(8);
        query.extend_from_slice(&genome[1000..1078]);
        query.extend_from_slice(&genome[3000..3064]);
        query[10] = b'C';
        for strand in [Strand::Plus, Strand::Minus] {
            let junction = Junction {
                contig: 0,
                start: 1078,
                end: 3000,
                minus_strand: true,
            };
            let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap())
                .with_lib(Arc::new(JunctionLib::build([junction])));
            let read = if strand == Strand::Plus {
                query.clone()
            } else {
                revcomp(&query)
            };
            let anchors = [(11, 1003, 72), (86, 3000, 64)]
                .into_iter()
                .flat_map(|(q, r, span)| {
                    (0..span - 15)
                        .step_by(4)
                        .map(move |offset| crate::seed::Anchor {
                            contig: 0,
                            strand,
                            qpos: q + offset,
                            rpos: r + offset,
                        })
                })
                .collect();
            let chain = Chain {
                contig: 0,
                strand,
                anchors,
                score: 90,
            };
            let hit = aligner.try_spliced(&read, &revcomp(&read), &chain).unwrap();
            assert_eq!(hit.pos, 1000);
            assert_eq!(
                hit.cigar,
                vec![
                    CigarOp::SoftClip(8),
                    CigarOp::Match(78),
                    CigarOp::RefSkip(1922),
                    CigarOp::Match(64),
                ]
            );
            let mut full = Placement::from_genomic(
                0,
                hit.pos,
                strand,
                hit.cigar,
                hit.junctions,
                90,
                hit.score,
            );
            let mut clipped = Placement::from_genomic(
                0,
                1003,
                strand,
                vec![
                    CigarOp::SoftClip(11),
                    CigarOp::Match(75),
                    CigarOp::RefSkip(1922),
                    CigarOp::Match(64),
                ],
                full.junctions.clone(),
                90,
                0,
            );
            aligner.rescore_strict(&query, &mut full);
            aligner.rescore_strict(&query, &mut clipped);
            assert_eq!(full.strict.score, clipped.strict.score);
            assert_eq!(full.strict.ea_tc, 1);
            assert_eq!(clipped.strict.ea_tc, 0);
            for mut candidates in [
                vec![full.clone(), clipped.clone()],
                vec![clipped.clone(), full.clone()],
            ] {
                dedup_by_footprint(&mut candidates, read.len());
                let result = aligner.arbitrate(&read, candidates, 10, true).unwrap();
                assert_eq!(result.pos, full.pos);
                assert_eq!(result.cigar, full.cigar);
                assert_eq!(result.ea_tc, 1);
            }
        }
    }

    #[test]
    fn intron_pool_proposal_survives_three_segment_best_chain() {
        // The intron constructor searches the complete anchor pool. Repeat
        // anchors in the best chain must not disable an independent exon pair.
        for strand in [Strand::Plus, Strand::Minus] {
            for leading_clip in [false, true] {
                let mut state = 0x719ac57bu64;
                let mut genome: Vec<u8> = (0..6000)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        b"ACGT"[(state & 3) as usize]
                    })
                    .collect();
                genome[2064..2066].copy_from_slice(b"GT");
                genome[3998..4000].copy_from_slice(b"AG");
                genome[2063] = b'T';
                genome[4000] = b'A';
                let mut query = Vec::new();
                if leading_clip {
                    query.extend_from_slice(&b"N".repeat(15));
                }
                query.extend_from_slice(&genome[2000..2064]);
                query.extend_from_slice(&genome[4000..4071]);
                if !leading_clip {
                    query.extend_from_slice(&b"N".repeat(15));
                }
                let offset = if leading_clip { 15 } else { 0 };
                genome[100..115].copy_from_slice(&query[offset + 5..offset + 20]);
                genome[1000..1030].copy_from_slice(&query[offset + 30..offset + 60]);
                let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                let read = if strand == Strand::Plus {
                    query.clone()
                } else {
                    revcomp(&query)
                };
                let chain = Chain {
                    contig: 0,
                    strand,
                    anchors: [
                        (offset + 5, 100),
                        (offset + 30, 1000),
                        (offset + 45, 1015),
                        (offset + 70, 4006),
                    ]
                    .map(|(q, r)| crate::seed::Anchor {
                        contig: 0,
                        strand,
                        qpos: q as u32,
                        rpos: r,
                    })
                    .to_vec(),
                    score: 62,
                };
                assert_eq!(segment_chain(&chain, 15, 20).len(), 3);
                let mins = crate::seed::minimizers(&read, aligner.config.seed);
                let mut all = Vec::new();
                collect_anchors_into(aligner.index, &mins, 150, 15, MID_OCC_CAP, &mut all);
                all.retain(|a| a.strand == strand);
                let hit = crate::intron_chain::try_intron_chain_placement(
                    aligner.index.reference,
                    &JunctionLib::default(),
                    &query,
                    &all,
                    15,
                    &crate::intron_chain::IntronParams::default(),
                    &ExtendParams::default(),
                    &mut ExtendBuffer::default(),
                )
                .expect("complete pool independently supports the exon pair");
                assert_eq!(hit.genomic_start(), 2000);
                let expected = if leading_clip {
                    vec![
                        CigarOp::SoftClip(15),
                        CigarOp::Match(64),
                        CigarOp::RefSkip(1936),
                        CigarOp::Match(71),
                    ]
                } else {
                    vec![
                        CigarOp::Match(64),
                        CigarOp::RefSkip(1936),
                        CigarOp::Match(71),
                        CigarOp::SoftClip(15),
                    ]
                };
                assert_eq!(hit.cigar, expected);
                let mut candidates = Vec::new();
                let mut budget = EXTEND_BUDGET;
                aligner.extend_chain(
                    &read,
                    &mins,
                    &chain,
                    true,
                    false,
                    &revcomp(&read),
                    &mut budget,
                    false,
                    &mut candidates,
                );
                assert!(
                    candidates
                        .iter()
                        .any(|p| p.pos == 2000 && p.cigar == expected),
                    "three-segment best chain suppressed the independent pool proposal"
                );
                assert!(
                    candidates.iter().any(|p| p.junctions.is_empty()),
                    "pool proposal must preserve continuous competing explanations"
                );
                for p in &mut candidates {
                    aligner.rescore_strict(&query, p);
                }
                dedup_by_footprint(&mut candidates, read.len());
                sort_placements(&mut candidates);
                let selected = aligner
                    .arbitrate(&read, candidates, mins.len(), true)
                    .expect("supported 135bp body clears the unchanged 99bp gate");
                assert_eq!(selected.pos, 2000);
                assert_eq!(selected.cigar, expected);
            }
        }
    }

    #[test]
    fn intron_fast_path_reconstructs_asymmetric_exons_with_technical_clips() {
        let mut state = 0x517cc1b7u64;
        let mut genome: Vec<u8> = (0..5000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        for (donor, acceptor) in [(b"GT", b"AG"), (b"CT", b"AC")] {
            genome[1060..1062].copy_from_slice(donor);
            genome[2998..3000].copy_from_slice(acceptor);
            // Unique seam in either motif orientation (no one-base slide).
            genome[1059] = b'T';
            genome[3000] = b'A';
            let mut query = b"N".repeat(11);
            query.extend_from_slice(&genome[1000..1060]);
            query.extend_from_slice(&genome[3000..3079]);
            let reference = ref_from(std::str::from_utf8(&genome).unwrap());
            for strand in [Strand::Plus, Strand::Minus] {
                let anchors: Vec<_> = [(11, 1000, 60), (71, 3000, 79)]
                    .into_iter()
                    .flat_map(|(q, r, span)| {
                        (3..span - 14)
                            .step_by(4)
                            .map(move |offset| crate::seed::Anchor {
                                contig: 0,
                                strand,
                                qpos: q + offset,
                                rpos: r + offset,
                            })
                    })
                    .collect();
                let hit = crate::intron_chain::try_intron_chain_placement(
                    &reference,
                    &JunctionLib::default(),
                    &query,
                    &anchors,
                    15,
                    &crate::intron_chain::IntronParams::default(),
                    &ExtendParams::default(),
                    &mut ExtendBuffer::default(),
                )
                .expect("supported fragment pair must retain an independent fast proposal");
                let p = Placement::from_genomic(
                    0,
                    hit.genomic_start(),
                    strand,
                    hit.cigar,
                    hit.junctions,
                    90,
                    hit.extension.score,
                );
                assert!(
                    p.has_consistent_junctions(),
                    "fast proposal CIGAR and junction disagree"
                );
                assert_eq!(p.pos, 1000);
                assert_eq!(
                    p.cigar,
                    vec![
                        CigarOp::SoftClip(11),
                        CigarOp::Match(60),
                        CigarOp::RefSkip(1940),
                        CigarOp::Match(79),
                    ]
                );
                assert_eq!(hit.extension.score, 278);
            }
        }
    }

    #[test]
    fn intron_proposal_junctions_must_match_cigar_reference_intervals() {
        let mut p = placement(
            0,
            100,
            &[
                CigarOp::SoftClip(5),
                CigarOp::Match(20),
                CigarOp::Ins(3),
                CigarOp::Del(2),
                CigarOp::RefSkip(50),
                CigarOp::Match(25),
            ],
        );
        p.junctions.push(RefinedJunction {
            junction: Junction {
                contig: 0,
                start: 122,
                end: 172,
                minus_strand: false,
            },
            known_support: 1,
            signal: SpliceSignal::GtAg,
        });
        assert!(p.has_consistent_junctions());
        p.junctions[0].junction.start += 1;
        assert!(!p.has_consistent_junctions());
        p.junctions[0].junction.start -= 1;
        p.junctions[0].junction.contig = 1;
        assert!(!p.has_consistent_junctions());
        p.junctions.clear();
        assert!(!p.has_consistent_junctions());
        assert!(placement(0, 100, &[CigarOp::Match(25)]).has_consistent_junctions());
    }

    #[test]
    fn intron_fast_path_keeps_generic_splicing_and_consistent_seams() {
        let mut state = 0x243f6a88u64;
        let mut random = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                })
                .collect()
        };
        let mut genome = random(8000);
        let exons = [(1000, 1032), (4257, 4333), (6144, 6175)];
        genome[1032..1034].copy_from_slice(b"GT");
        genome[4255..4257].copy_from_slice(b"AG");
        genome[4333..4335].copy_from_slice(b"GT");
        genome[6142..6144].copy_from_slice(b"AG");
        let mut query = b"N".repeat(11);
        for (start, end) in exons {
            query.extend_from_slice(&genome[start..end]);
        }
        assert_eq!(query.len(), 150);
        for strand in [Strand::Plus, Strand::Minus] {
            let library = JunctionLib::build([
                Junction {
                    contig: 0,
                    start: 1032,
                    end: 4257,
                    minus_strand: strand == Strand::Minus,
                },
                Junction {
                    contig: 0,
                    start: 4333,
                    end: 6144,
                    minus_strand: strand == Strand::Minus,
                },
            ]);
            let mut aligner =
                aligner_on(std::str::from_utf8(&genome).unwrap()).with_lib(Arc::new(library));
            let read = if strand == Strand::Plus {
                query.clone()
            } else {
                revcomp(&query)
            };
            let mut anchors = Vec::new();
            for (q, r, span) in [(43, 4257, 76), (119, 6144, 31)] {
                for offset in (3..span - 14).step_by(4) {
                    anchors.push(crate::seed::Anchor {
                        contig: 0,
                        strand,
                        qpos: q + offset,
                        rpos: r + offset,
                    });
                }
            }
            let chain = Chain {
                contig: 0,
                strand,
                anchors,
                score: 93,
            };
            let expected = aligner
                .try_spliced(&read, &revcomp(&read), &chain)
                .expect("generic path must have an independent supported explanation");
            assert_eq!(expected.pos, 4257);
            assert_eq!(
                expected.cigar,
                vec![
                    CigarOp::SoftClip(43),
                    CigarOp::Match(76),
                    CigarOp::RefSkip(1811),
                    CigarOp::Match(31),
                ]
            );
            let mut candidates = Vec::new();
            aligner.extend_chain(
                &read,
                &crate::seed::minimizers(&read, aligner.config.seed),
                &chain,
                true,
                false,
                &revcomp(&read),
                &mut 8,
                false,
                &mut candidates,
            );
            assert!(
                candidates
                    .iter()
                    .any(|p| p.pos == expected.pos && p.cigar == expected.cigar),
                "an intron fast-path proposal must not suppress the generic splice candidate"
            );
            for p in &candidates {
                let mut rf = p.pos;
                for op in &p.cigar {
                    match op {
                        CigarOp::Match(n) | CigarOp::Del(n) => rf += *n,
                        CigarOp::RefSkip(n) => {
                            assert!(
                                p.junctions
                                    .iter()
                                    .any(|j| j.junction.start == rf && j.junction.end == rf + n),
                                "published CIGAR seam must agree with its refined junction"
                            );
                            rf += *n;
                        }
                        _ => {}
                    }
                }
            }
        }
    }

    #[test]
    fn hyper_short_trigger_rescores_every_competing_locus_in_one_domain() {
        let mut genome = b"C".repeat(3000);
        let mut query = b"A".repeat(150);
        for q in [5, 15, 25, 35, 45, 55] {
            query[q] = b'C';
        }
        for (q, base) in query.iter().enumerate() {
            genome[100 + q] = *base;
            genome[1100 + q] = *base;
        }
        // Each full candidate has six T>C edits (4%, below trigger density).
        for q in [5, 15, 25, 35, 45, 55] {
            genome[100 + q] = b'T';
            genome[1100 + q] = b'T';
        }
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut candidates = Vec::new();
        for pos in [100, 1100] {
            let mut p = placement(0, pos, &[CigarOp::Match(150)]);
            p.origin = crate::mapq::Origin::L1; // no reextension, isolates common-domain scoring
            aligner.rescore_strict(&query, &mut p);
            assert!(hyper_direction(&p).is_none());
            candidates.push(p);
        }
        let mut short = placement(0, 100, &[CigarOp::Match(84), CigarOp::SoftClip(66)]);
        short.origin = crate::mapq::Origin::L1;
        aligner.rescore_strict(&query, &mut short);
        assert_eq!(hyper_direction(&short), Some(EditDir::Tc));
        candidates.push(short);
        dedup_by_footprint(&mut candidates, query.len());
        aligner.hyper_rearbitrate(&query, &revcomp(&query), &mut candidates);
        assert!(candidates.iter().all(|p| p.ea_dir == EditDir::Tc));
        let result = aligner.arbitrate(&query, candidates, 10, true).unwrap();
        assert_eq!(result.score, 288);
        assert_eq!(result.second_score, 288);
        assert_eq!(result.mapq_of_checked(), 0);
        assert_eq!(result.n_cands, 2);
    }

    fn hyper_opposite_strand_case(trigger_strand: Strand, trigger_dir: EditDir, genomic: bool) {
        // Purely synthetic two-locus competition. The strict leader has seven
        // edits (below the 5% trigger), while its opposite-strand competitor
        // has nine edits plus an ordinary mismatch and triggers re-arbitration.
        // A read-oriented conversion must complement when the candidate flips.
        let mut state = 0x91cb_ab67u64;
        let mut random = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                })
                .collect()
        };
        let mut read = random(150);
        let edited_positions = [15, 27, 39, 51, 63, 75, 87, 99, 111];
        let trigger_base = if trigger_dir == EditDir::Ag {
            b'G'
        } else {
            b'C'
        };
        let raw_base = if trigger_strand == Strand::Plus {
            trigger_base
        } else if trigger_base == b'G' {
            b'C'
        } else {
            b'G'
        };
        for &i in &edited_positions {
            read[i] = raw_base;
        }
        let other_strand = if trigger_strand == Strand::Plus {
            Strand::Minus
        } else {
            Strand::Plus
        };
        let other_dir = if trigger_dir == EditDir::Ag {
            EditDir::Tc
        } else {
            EditDir::Ag
        };
        let rc = revcomp(&read);
        let oriented = |strand| {
            if strand == Strand::Plus {
                read.clone()
            } else {
                rc.clone()
            }
        };
        let mut trigger_reference = oriented(trigger_strand);
        let mut other_reference = oriented(other_strand);
        for (number, &i) in edited_positions.iter().enumerate() {
            let t = if trigger_strand == Strand::Plus {
                i
            } else {
                149 - i
            };
            let o = if other_strand == Strand::Plus {
                i
            } else {
                149 - i
            };
            trigger_reference[t] = if trigger_dir == EditDir::Ag {
                b'A'
            } else {
                b'T'
            };
            if number < 7 {
                other_reference[o] = if other_dir == EditDir::Ag { b'A' } else { b'T' };
            }
        }
        let ordinary = if trigger_strand == Strand::Plus {
            130
        } else {
            19
        };
        let qb = trigger_reference[ordinary];
        trigger_reference[ordinary] = *b"ACGT"
            .iter()
            .find(|&&rb| rb != qb && crate::extend::exempt_dir(rb, qb).is_none())
            .unwrap();
        let mut genome = random(1200);
        genome[100..250].copy_from_slice(&trigger_reference);
        genome[800..950].copy_from_slice(&other_reference);
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut trigger = placement(0, 100, &[CigarOp::Match(150)]);
        trigger.strand = trigger_strand;
        let mut other = placement(0, 800, &[CigarOp::Match(150)]);
        other.strand = other_strand;
        for p in [&mut trigger, &mut other] {
            if !genomic {
                p.origin = crate::mapq::Origin::L1;
            }
            let query = if p.strand == Strand::Plus { &read } else { &rc };
            aligner.rescore_strict(query, p);
        }
        assert_eq!(trigger.strict.score, 240);
        assert_eq!(other.strict.score, 258);
        assert_eq!(hyper_direction(&trigger), Some(trigger_dir));
        assert_eq!(hyper_direction(&other), None);
        let mut candidates = vec![other, trigger];
        sort_placements(&mut candidates);
        aligner.hyper_rearbitrate(&read, &rc, &mut candidates);
        let other = candidates.iter().find(|p| p.pos == 800).unwrap();
        assert_eq!(
            other.ea_dir, other_dir,
            "opposite strand must complement the read's domain"
        );
        assert_eq!(other.strict.score, 286);
        let trigger = candidates.iter().find(|p| p.pos == 100).unwrap();
        assert_eq!(trigger.ea_dir, trigger_dir);
        assert_eq!(trigger.strict.score, 276);
        let result = aligner.arbitrate(&read, candidates, 10, true).unwrap();
        assert_eq!(result.pos, 800);
        assert_eq!(result.strand, other_strand);
        assert_eq!(result.second_score, 276);
    }

    #[test]
    fn hyper_opposite_strand_candidates_share_one_read_oriented_domain() {
        for strand in [Strand::Plus, Strand::Minus] {
            for dir in [EditDir::Ag, EditDir::Tc] {
                hyper_opposite_strand_case(strand, dir, false);
            }
        }
    }

    #[test]
    fn hyper_opposite_strand_genomic_reextension_uses_complemented_domain() {
        for strand in [Strand::Plus, Strand::Minus] {
            for dir in [EditDir::Ag, EditDir::Tc] {
                hyper_opposite_strand_case(strand, dir, true);
            }
        }
    }

    #[test]
    fn hyper_complemented_reextension_recovers_an_edited_clipped_flank() {
        // The opposite-strand trigger fixes the read's conversion alphabet.
        // Its genomic competitor starts with an unedited 129 bp prefix;
        // recovering the last 21 bp requires exempting nine edits there.
        // Common rescoring alone cannot restore those missing CIGAR columns.
        for other_strand in [Strand::Plus, Strand::Minus] {
            for other_dir in [EditDir::Ag, EditDir::Tc] {
                let mut state = 0x98b3_5f2du64;
                let mut random = |len: usize| -> Vec<u8> {
                    (0..len)
                        .map(|_| {
                            state ^= state << 13;
                            state ^= state >> 7;
                            state ^= state << 17;
                            b"ACGT"[(state & 3) as usize]
                        })
                        .collect()
                };
                let mut query = random(150);
                let edited_base = if other_dir == EditDir::Ag { b'G' } else { b'C' };
                let reference_base = if other_dir == EditDir::Ag { b'A' } else { b'T' };
                for i in 0..9 {
                    query[129 + 2 * i] = edited_base;
                }
                let mut other_reference = query.clone();
                for i in 0..9 {
                    other_reference[129 + 2 * i] = reference_base;
                }
                let trigger_query = revcomp(&query);
                let mut trigger_reference = trigger_query.clone();
                for i in 0..9 {
                    trigger_reference[149 - (129 + 2 * i)] =
                        if other_dir == EditDir::Ag { b'T' } else { b'A' };
                }
                let qb = trigger_reference[100];
                trigger_reference[100] = *b"ACGT"
                    .iter()
                    .find(|&&rb| rb != qb && crate::extend::exempt_dir(rb, qb).is_none())
                    .unwrap();
                let mut genome = random(1200);
                genome[100..250].copy_from_slice(&trigger_reference);
                genome[800..950].copy_from_slice(&other_reference);
                let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
                let read = if other_strand == Strand::Plus {
                    query.clone()
                } else {
                    revcomp(&query)
                };
                let mut trigger = placement(0, 100, &[CigarOp::Match(150)]);
                trigger.strand = if other_strand == Strand::Plus {
                    Strand::Minus
                } else {
                    Strand::Plus
                };
                trigger.origin = crate::mapq::Origin::L1;
                aligner.rescore_strict(&trigger_query, &mut trigger);
                let mut other = placement(0, 800, &[CigarOp::Match(129), CigarOp::SoftClip(21)]);
                other.strand = other_strand;
                aligner.rescore_strict(&query, &mut other);
                assert_eq!(other.strict.score, 258);
                assert_eq!(hyper_direction(&other), None);
                let mut candidates = vec![other, trigger];
                aligner.hyper_rearbitrate(&read, &revcomp(&read), &mut candidates);
                let recovered = candidates.iter().find(|p| p.pos == 800).unwrap();
                assert_eq!(recovered.cigar, vec![CigarOp::Match(150)]);
                assert_eq!(recovered.ea_dir, other_dir);
                assert_eq!(recovered.strict.score, 282);
                assert_eq!(candidates[0].pos, 800);
            }
        }
    }

    #[test]
    fn hyper_same_strand_competitor_cannot_choose_its_own_editing_direction() {
        let read = b"CG".repeat(75);
        let mut trigger_ref = read.clone();
        let mut other_ref = read.clone();
        for i in 0..9 {
            trigger_ref[21 + 10 * i] = b'A';
        }
        trigger_ref[130] = b'G'; // ordinary ref G/read C
        for i in 0..7 {
            other_ref[20 + 10 * i] = b'T';
        }
        let mut genome = b"N".repeat(1200);
        genome[100..250].copy_from_slice(&trigger_ref);
        genome[800..950].copy_from_slice(&other_ref);
        let mut aligner = aligner_on(std::str::from_utf8(&genome).unwrap());
        let mut trigger = placement(0, 100, &[CigarOp::Match(150)]);
        let mut other = placement(0, 800, &[CigarOp::Match(150)]);
        for p in [&mut trigger, &mut other] {
            p.origin = crate::mapq::Origin::L1;
            aligner.rescore_strict(&read, p);
        }
        let mut candidates = vec![other, trigger];
        sort_placements(&mut candidates);
        aligner.hyper_rearbitrate(&read, &revcomp(&read), &mut candidates);
        assert!(candidates.iter().all(|p| p.ea_dir == EditDir::Ag));
        assert_eq!(
            candidates
                .iter()
                .find(|p| p.pos == 800)
                .unwrap()
                .strict
                .score,
            258
        );
        assert_eq!(candidates[0].pos, 100);
    }

    #[test]
    fn refuse_placement_below_two_thirds_coverage() {
        // 300bp contig; read = 75bp of the contig + 75bp absent sequence.
        let contig = "GATTACAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGCAAGTGGCATTGACCAAGTT"
            .to_string()
            + "TGCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATT";
        let mut aln = aligner_on(&contig);
        let head = &contig.as_bytes()[..75];
        let tail = b"T".repeat(75).leak();
        let mut read = head.to_vec();
        read.extend_from_slice(&tail[..]);
        assert_eq!(read.len(), 150);
        assert!(aln.align_read(&read).is_none());
    }

    #[test]
    fn repeat_copy_scores_mapq_zero() {
        // Two identical 120bp copies: a 60bp read from the copy seeds both.
        let copy = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGACCAAGTTACAGGATCCGATTAAGGC";
        let spacer = "A".repeat(120);
        let seq = format!("{copy}{spacer}{copy}");
        let mut aln = aligner_on(&seq);
        let read = &copy.as_bytes()[30..90];
        let aln_out = aln.align_read(read).expect("placed");
        assert_eq!(aln_out.mapq_of_checked(), 0);
    }

    /// Exon1 | GT..AG intron(60) | exon2 fixture with a 40bp read spanning
    /// 20+20 exon bases and one 15-mer anchor per exon (read-gap 5,
    /// ref-gap 65 ⇒ implied intron exactly 60).
    fn spliced_chain_fixture() -> (String, Vec<u8>, Chain) {
        let exon1 = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCAAGTTGACGTTACAGGCATTAC"; // 57
        let intron = format!("GT{}AG", "N".repeat(56)); // 60, GT..AG inside
        let exon2 = "CATGGACTTGAACCTTGACGTACATGGACTTGAACCTTGACGTA"; // 44
        let seq = format!("{exon1}{intron}{exon2}");
        let s = seq.as_bytes();
        let start = (exon1.len() - 20) as u32; // last 20 exon1 bases
        let e2 = exon1.len() + intron.len();
        let mut read: Vec<u8> = s[start as usize..exon1.len()].to_vec();
        read.extend_from_slice(&s[e2..e2 + 20]);
        let chain = Chain {
            contig: 0,
            strand: Strand::Plus,
            anchors: vec![
                crate::seed::Anchor {
                    qpos: 0,
                    rpos: start,
                    contig: 0,
                    strand: Strand::Plus,
                },
                crate::seed::Anchor {
                    qpos: 20,
                    rpos: e2 as u32,
                    contig: 0,
                    strand: Strand::Plus,
                },
            ],
            score: 200,
        };
        (seq, read, chain)
    }

    #[test]
    fn try_spliced_stitches_two_segment_chain() {
        let (seq, read, chain) = spliced_chain_fixture();
        let mut aln = aligner_on(&seq);
        let rc = revcomp(&read);
        let hit = aln.try_spliced(&read, &rc, &chain).expect("spliced hit");
        assert_eq!(hit.pos, 37); // exon1 is 57bp: last 20 exon1 bases start at 37
        assert_eq!(
            hit.cigar,
            vec![CigarOp::Match(20), CigarOp::RefSkip(60), CigarOp::Match(20)]
        );
        assert_eq!(hit.junctions.len(), 1);
        assert_eq!(hit.junctions[0].signal, SpliceSignal::GtAg);
        assert_eq!(hit.junctions[0].junction.start, 57);
        assert_eq!(hit.junctions[0].junction.end, 117);
    }
}

#[cfg(test)]
pub(crate) mod background_tests;
