//! Observable qualification of coverage-ineligible paired RNA candidates.
//! Indexed sketch support is not exhaustive uniqueness; output is always MAPQ0.
use super::*;
use crate::diagnostics::{
    BoundaryCoreProposal, ShortBodyCandidate, ShortBodyDecision, ShortBodySketchAttempt,
};
use crate::index::unpack_pos;

/// Exact, complex body with one observably separated canonical terminal tail.
/// Query is already in candidate alignment orientation.
fn exact_body_seeds(
    index: &Index,
    query: &[u8],
    p: &Placement,
) -> Option<Vec<(u32, crate::seed::Minimizer)>> {
    if p.primary_eligible(query.len()) || p.strict.score <= 0 || !p.has_consistent_junctions() {
        return None;
    }
    let reference = index.reference.contigs.get(p.contig as usize)?;
    let mut clip = None;
    for (i, op) in p.cigar.iter().enumerate() {
        match op {
            CigarOp::SoftClip(n) if (i == 0 || i + 1 == p.cigar.len()) && clip.is_none() => {
                clip = Some((i == 0, *n as usize));
            }
            CigarOp::Match(n) | CigarOp::RefSkip(n) if *n > 0 => {}
            _ => return None,
        }
    }
    let (prefix, n) = clip?;
    if n < 30 || n >= query.len() {
        return None;
    }
    let (tail, adjacent) = if prefix {
        (&query[..n], query[n])
    } else {
        (&query[query.len() - n..], query[query.len() - n - 1])
    };
    let base = tail[0].to_ascii_uppercase();
    if !b"ACGT".contains(&base)
        || !tail.iter().all(|b| b.to_ascii_uppercase() == base)
        || !b"ACGT".contains(&adjacent.to_ascii_uppercase())
        || adjacent.to_ascii_uppercase() == base
        || !matches!(
            p.cigar.first(),
            Some(CigarOp::Match(_) | CigarOp::SoftClip(_))
        )
        || !matches!(
            p.cigar.last(),
            Some(CigarOp::Match(_) | CigarOp::SoftClip(_))
        )
    {
        return None;
    }
    let mut q = 0usize;
    let mut r = p.pos;
    let mut counts = [0usize; 4];
    let mut seeds = Vec::new();
    for op in &p.cigar {
        match *op {
            CigarOp::Match(n) => {
                let qe = q.checked_add(n as usize)?;
                let re = r.checked_add(n)?;
                let block = query.get(q..qe)?;
                if re > reference.len {
                    return None;
                }
                for (i, b) in block.iter().enumerate() {
                    let b = b.to_ascii_uppercase();
                    let code = b"ACGT".iter().position(|&a| a == b)?;
                    if reference.base(r + i as u32).to_ascii() != b {
                        return None;
                    }
                    counts[code] += 1;
                }
                for seed in crate::seed::minimizers(block, index.params) {
                    seeds.push((
                        q as u32 + seed.pos,
                        crate::seed::Minimizer {
                            pos: r + seed.pos,
                            ..seed
                        },
                    ));
                }
                q = qe;
                r = re;
            }
            CigarOp::SoftClip(n) => q = q.checked_add(n as usize)?,
            CigarOp::RefSkip(n) => {
                r = r.checked_add(n)?;
                if r > reference.len {
                    return None;
                }
            }
            _ => return None,
        }
    }
    let body: usize = counts.iter().sum();
    if q != query.len()
        || body < 31
        || counts.iter().filter(|&&n| n > 0).count() < 3
        || counts.iter().any(|&n| n * 5 > body * 4)
    {
        return None;
    }
    Some(seeds)
}

fn disjoint_support(mut supported: Vec<(u32, u64)>, k: u32) -> Vec<(u32, u64)> {
    supported.sort_unstable();
    supported.dedup();
    let mut end = 0;
    let mut disjoint = Vec::new();
    for (start, key) in supported {
        if start >= end {
            disjoint.push((start, key));
            end = start + k;
        }
    }
    disjoint
}

/// Original global-index support contract, unchanged by conditional recovery.
fn indexed_support(index: &Index, query: &[u8], p: &Placement) -> Option<Vec<(u32, u64)>> {
    let mut supported = Vec::new();
    for (query_position, seed) in exact_body_seeds(index, query, p)? {
        let hits = index.query(seed.kmer);
        if hits.count != 1 {
            continue;
        }
        let (contig, position, strand) = unpack_pos(hits.positions[0]);
        if contig == p.contig && position == seed.pos && strand == seed.strand {
            supported.push((query_position, seed.kmer));
        }
    }
    let disjoint = disjoint_support(supported, index.params.k);
    (disjoint.len() >= 2).then_some(disjoint)
}

/// Bounded supplied-index evidence conditioned on the actual complete mate.
/// Only contiguous bodies use this proposal; the original route handles N.
fn mate_conditioned_support(
    index: &Index,
    query: &[u8],
    p: &Placement,
    mate: &ReadAlignment,
    attempt: &mut ShortBodySketchAttempt,
) -> Option<Vec<(u32, u64)>> {
    if p.cigar
        .iter()
        .any(|op| !matches!(op, CigarOp::Match(_) | CigarOp::SoftClip(_)))
    {
        return None;
    }
    let mut seeds = exact_body_seeds(index, query, p)?;
    seeds.sort_unstable_by_key(|(q, seed)| (index.query(seed.kmer).count, *q, seed.kmer));
    let body = p
        .cigar
        .iter()
        .filter_map(|op| match op {
            CigarOp::Match(n) => Some(*n),
            _ => None,
        })
        .sum::<u32>();
    let mut supported = Vec::new();
    for (q, seed) in seeds {
        let hits = index.query(seed.kmer);
        if hits.count == 0 {
            continue;
        }
        if hits.count > MID_OCC_CAP {
            attempt.over_cap_seed_lists += 1;
            attempt.omitted_seed_lists += 1;
            continue;
        }
        if hits.count as usize > MID_ROUND_MAX_ANCHORS - attempt.scanned_postings {
            attempt.omitted_seed_lists += 1;
            continue;
        }
        attempt.scanned_seed_lists += 1;
        attempt.scanned_postings += hits.count as usize;
        let body_offset = seed.pos.checked_sub(p.pos)?;
        let mut conditioned_positions = 0usize;
        let mut expected_position = false;
        for &packed in hits.positions {
            let (contig, position, strand) = unpack_pos(packed);
            if contig != mate.contig || strand != seed.strand {
                continue;
            }
            let Some(start) = position.checked_sub(body_offset) else {
                continue;
            };
            let Some(end) = start.checked_add(body) else {
                continue;
            };
            if end > index.reference.contigs[contig as usize].len {
                continue;
            }
            let hypothetical = ReadAlignment {
                contig,
                pos: start,
                strand: p.strand,
                cigar: vec![CigarOp::Match(body)],
                ..Default::default()
            };
            if crate::pair::relate(&hypothetical, mate) == crate::pair::PairRelation::Proper {
                conditioned_positions += 1;
                expected_position |= start == p.pos;
            }
        }
        if conditioned_positions == 1 && expected_position {
            supported.push((q, seed.kmer));
        }
    }
    let disjoint = disjoint_support(supported, index.params.k);
    attempt.retained_disjoint_seeds = disjoint.len();
    (disjoint.len() >= 2).then_some(disjoint)
}

fn same_path(a: &Placement, b: &Placement) -> bool {
    a.contig == b.contig && a.pos == b.pos && a.strand == b.strand && a.cigar == b.cigar
}

/// Project only columns invariant to biological/technical boundary assignment.
/// Neither the run nor its soft clip is asserted to be technical sequence.
fn boundary_core(index: &Index, query: &[u8], p: &Placement) -> Option<(Placement, u32, u32, u8)> {
    if p.primary_eligible(query.len()) || p.strict.score <= 0 {
        return None;
    }
    let (prefix, clipped, aligned) = match p.cigar.as_slice() {
        [CigarOp::SoftClip(s), CigarOp::Match(m)] => (true, *s, *m),
        [CigarOp::Match(m), CigarOp::SoftClip(s)] => (false, *s, *m),
        _ => return None,
    };
    if clipped < 30 || aligned == 0 || clipped as usize + aligned as usize != query.len() {
        return None;
    }
    let clipped = clipped as usize;
    let aligned = aligned as usize;
    let (tail, body) = if prefix {
        (&query[..clipped], &query[clipped..])
    } else {
        (&query[aligned..], &query[..aligned])
    };
    let base = tail[0].to_ascii_uppercase();
    if !b"ACGT".contains(&base) || tail.iter().any(|b| b.to_ascii_uppercase() != base) {
        return None;
    }
    let reference = index.reference.contigs.get(p.contig as usize)?;
    let end = p.pos.checked_add(aligned as u32)?;
    if end > reference.len
        || body.iter().enumerate().any(|(i, &b)| {
            !b"ACGT".contains(&b.to_ascii_uppercase())
                || reference.base(p.pos + i as u32).to_ascii() != b.to_ascii_uppercase()
        })
    {
        return None;
    }
    let run = if prefix {
        body.iter()
            .take_while(|b| b.to_ascii_uppercase() == base)
            .count()
    } else {
        body.iter()
            .rev()
            .take_while(|b| b.to_ascii_uppercase() == base)
            .count()
    };
    if run == 0 || run == aligned {
        return None;
    }
    let mut core = p.clone();
    let n = (aligned - run) as u32;
    let s = (clipped + run) as u32;
    core.cigar = if prefix {
        core.pos = core.pos.checked_add(run as u32)?;
        vec![CigarOp::SoftClip(s), CigarOp::Match(n)]
    } else {
        vec![CigarOp::Match(n), CigarOp::SoftClip(s)]
    };
    Some((
        core,
        if prefix { clipped as u32 } else { n },
        run as u32,
        base,
    ))
}

impl Aligner<'_> {
    /// Pipeline-only fallback; callers must preserve every already mapped read.
    /// Full mate geometry is essential; coordinate-only legacy APIs do not use it.
    pub(crate) fn recover_short_body_observed(
        &mut self,
        seq: &[u8],
        mate: &ReadAlignment,
        cached: Option<&CachedPlacements>,
        diagnostic: Option<&mut Option<crate::diagnostics::CandidatePool>>,
    ) -> Option<ReadAlignment> {
        if !self.config.rna
            || !self.config.mate_short_body
            || mate.rescued
            || !(crate::mapq::PILEUP_GATE..255).contains(&crate::mapq::mapq_of(mate))
            || self
                .index
                .reference
                .contigs
                .get(mate.contig as usize)
                .is_none()
            || !mate
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::Match(n) if *n > 0))
        {
            return None;
        }
        let (mut placements, n_seeds, _, budget) = match cached {
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
        // Same bounded ordinary window as the old mate attempt. Keep all exact
        // paths here, rather than footprint-deduplicating before qualification.
        if let Some(win) = self.mate_window_placement(seq, mate.contig, mate.pos) {
            placements.push(win);
        }
        let rc = self.take_rc(seq);
        for p in &mut placements {
            let query = if p.strand == Strand::Minus { &rc } else { seq };
            self.rescore_strict(query, p);
        }
        let mut boundary_core_proposals = Vec::new();
        let original_count = placements.len();
        let mut qualifying: Vec<(usize, Vec<(u32, u64)>, &'static str)> = Vec::new();
        let mut conditional_sketch_attempts = Vec::new();
        for round in 0..2 {
            if round == 1 {
                // Preserve every previously accepted or ambiguously qualified
                // original decision. Partial-core evidence is a final fallback.
                if !self.config.boundary_core || !qualifying.is_empty() {
                    break;
                }
                for source_rank in 0..original_count {
                    let p = &placements[source_rank];
                    let query = if p.strand == Strand::Minus { &rc } else { seq };
                    if let Some((mut core, start, length, base)) =
                        boundary_core(self.index, query, p)
                    {
                        self.rescore_strict(query, &mut core);
                        boundary_core_proposals.push(BoundaryCoreProposal {
                            source_rank,
                            rank_in_pool: placements.len(),
                            uncertainty_query_start: start,
                            uncertainty_bases: length,
                            run_base: match base {
                                b'A' => "A",
                                b'C' => "C",
                                b'G' => "G",
                                _ => "T",
                            },
                            assignment: "unknown_biological_or_technical",
                        });
                        placements.push(core);
                    }
                }
            }
            for (rank, p) in
                placements
                    .iter()
                    .enumerate()
                    .skip(if round == 0 { 0 } else { original_count })
            {
                if p.contig != mate.contig || p.strand == mate.strand {
                    continue;
                }
                let query = if p.strand == Strand::Minus { &rc } else { seq };
                let geometry = ReadAlignment {
                    contig: p.contig,
                    pos: p.pos,
                    strand: p.strand,
                    cigar: p.cigar.clone(),
                    ..Default::default()
                };
                if crate::pair::relate(&geometry, mate) != crate::pair::PairRelation::Proper {
                    continue;
                }
                let (seeds, scope) = if let Some(seeds) = indexed_support(self.index, query, p) {
                    (seeds, "global_index")
                } else {
                    let mut attempt = ShortBodySketchAttempt {
                        rank_in_pool: rank,
                        ..Default::default()
                    };
                    let seeds = mate_conditioned_support(self.index, query, p, mate, &mut attempt);
                    conditional_sketch_attempts.push(attempt);
                    let Some(seeds) = seeds else { continue };
                    (seeds, "mate_conditioned_index")
                };
                if !qualifying
                    .iter()
                    .any(|(i, _, _)| same_path(&placements[*i], p))
                {
                    qualifying.push((rank, seeds, scope));
                }
            }
        }
        if let Some(slot) = diagnostic {
            let mut pool = crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "mate_short_body",
            );
            pool.seed_budget = budget;
            pool.boundary_core_proposals = boundary_core_proposals.clone();
            pool.short_body = Some(ShortBodyDecision {
                query_coordinate_frame: "raw_alignment_query",
                decision: match qualifying.len() {
                    0 => "no_qualified_path",
                    1 => "qualified_mapq0",
                    _ => "ambiguous_paths",
                },
                qualified_paths: qualifying.len(),
                supported_candidates: qualifying
                    .iter()
                    .map(|(rank, seeds, scope)| ShortBodyCandidate {
                        rank_in_pool: *rank,
                        support_scope: scope,
                        disjoint_unique_indexed_seeds: if *scope == "global_index" {
                            seeds.clone()
                        } else {
                            Vec::new()
                        },
                        disjoint_mate_conditioned_indexed_seeds: if *scope
                            == "mate_conditioned_index"
                        {
                            seeds.clone()
                        } else {
                            Vec::new()
                        },
                    })
                    .collect(),
                conditional_sketch_attempts,
            });
            *slot = Some(pool);
        }
        self.scratch.rc = rc;
        if qualifying.len() != 1 {
            return None;
        }
        let (rank, seeds, _) = &qualifying[0];
        let p = &placements[*rank];
        Some(ReadAlignment {
            contig: p.contig,
            pos: p.pos,
            strand: p.strand,
            cigar: p.cigar.clone(),
            score: p.strict.score,
            second_score: placements
                .iter()
                .filter(|other| !same_path(other, p))
                .map(|other| other.strict.score)
                .filter(|&s| s > 0)
                .max()
                .unwrap_or(0),
            chain_score: p.chain_score,
            n_anchors: seeds.len(),
            junctions: p.junctions.clone(),
            n_seeds,
            rescued: true,
            short_body_qualified: true,
            pair_substitution_qualified: false,
            boundary_core_uncertainty: boundary_core_proposals
                .iter()
                .filter(|proposal| same_path(&placements[proposal.rank_in_pool], p))
                .fold(None, |range, proposal| {
                    let start = proposal.uncertainty_query_start;
                    let end = start + proposal.uncertainty_bases;
                    Some(match range {
                        None => (start, end - start),
                        Some((old_start, length)) => {
                            let left = start.min(old_start);
                            (left, end.max(old_start + length) - left)
                        }
                    })
                }),
            mapq_override: Some(0),
            origin: p.origin,
            n_cands: 1,
            seed_search_limited: budget.is_some_and(|b| b.limited()),
            ..Default::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index(genome: &[u8]) -> Index<'static> {
        Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), genome, b"\n"].concat()).unwrap(),
            SeedParams { k: 15, w: 5 },
        )
    }

    fn candidate(pos: u32, strand: Strand, cigar: Vec<CigarOp>) -> Placement {
        let mut p = Placement::from_genomic(0, pos, strand, cigar, Vec::new(), 30, 160);
        let mut r = pos;
        for op in &p.cigar {
            match *op {
                CigarOp::Match(n) => {
                    p.strict.aligned += n;
                    r += n;
                }
                CigarOp::RefSkip(n) => {
                    p.junctions.push(RefinedJunction {
                        junction: Junction {
                            contig: 0,
                            start: r,
                            end: r + n,
                            minus_strand: strand == Strand::Minus,
                        },
                        signal: SpliceSignal::NonCanonical,
                        known_support: 0,
                    });
                    r += n;
                }
                _ => {}
            }
        }
        p.strict.score = 2 * p.strict.aligned as i32;
        p
    }

    fn cache(placements: Vec<Placement>) -> CachedPlacements {
        CachedPlacements {
            placements,
            n_seeds: 12,
            genomic_any: true,
            seed_budget: None,
        }
    }

    fn mate(pos: u32, strand: Strand) -> ReadAlignment {
        ReadAlignment {
            contig: 0,
            pos,
            strand,
            cigar: vec![CigarOp::Match(150)],
            mapq_override: Some(20),
            ..Default::default()
        }
    }

    fn fixture() -> (Vec<u8>, Vec<u8>, Placement) {
        let genome = crate::fragment::tests::body(8000);
        let body = &genome[2000..2080];
        let tail = *b"ACGT"
            .iter()
            .find(|&&b| b != body[79] && b != genome[2080])
            .unwrap();
        let seq = [body, &vec![tail; 70]].concat();
        (
            genome,
            seq,
            candidate(
                2000,
                Strand::Plus,
                vec![CigarOp::Match(80), CigarOp::SoftClip(70)],
            ),
        )
    }

    #[test]
    fn boundary_core_identical_worlds_share_only_complex_columns() {
        for &base in b"ACGT" {
            for run in [1usize, 4, 7, 31] {
                for prefix in [false, true] {
                    let mut genome = crate::fragment::tests::body(8000);
                    let different = *b"ACGT".iter().find(|&&b| b != base).unwrap();
                    if prefix {
                        genome[2000..2000 + run].fill(base);
                        genome[2000 + run] = different;
                    } else {
                        genome[2080 - run..2080].fill(base);
                        genome[2079 - run] = different;
                    }
                    let body = genome[2000..2080].to_vec();
                    let oriented = if prefix {
                        [&vec![base; 70], body.as_slice()].concat()
                    } else {
                        [body.as_slice(), &vec![base; 70]].concat()
                    };
                    let idx = index(&genome);
                    for strand in [Strand::Plus, Strand::Minus] {
                        let p = candidate(
                            2000,
                            strand,
                            if prefix {
                                vec![CigarOp::SoftClip(70), CigarOp::Match(80)]
                            } else {
                                vec![CigarOp::Match(80), CigarOp::SoftClip(70)]
                            },
                        );
                        assert!(
                            indexed_support(&idx, &oriented, &p).is_none(),
                            "the biological/technical boundary remains unidentified"
                        );
                        let raw = if strand == Strand::Minus {
                            revcomp(&oriented)
                        } else {
                            oriented.clone()
                        };
                        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
                        let mut diagnostic = None;
                        let result = aligner.recover_short_body_observed(
                            &raw,&mate(2300,strand.flip()),Some(&cache(vec![p])),
                            Some(&mut diagnostic),
                        ).expect("invariant complex columns can emit without identifying the run boundary");
                        let expected = if prefix {
                            vec![
                                CigarOp::SoftClip((70 + run) as u32),
                                CigarOp::Match((80 - run) as u32),
                            ]
                        } else {
                            vec![
                                CigarOp::Match((80 - run) as u32),
                                CigarOp::SoftClip((70 + run) as u32),
                            ]
                        };
                        assert_eq!(result.cigar, expected);
                        assert_eq!(result.pos, 2000 + if prefix { run as u32 } else { 0 });
                        assert_eq!(crate::mapq::mapq_of(&result), 0);
                        // Every possible assignment of the same observed run to
                        // biological sequence contains the emitted core. These
                        // labels never enter recover_short_body_observed.
                        let core_start = if prefix { 70 + run } else { 0 };
                        let core_end = core_start + 80 - run;
                        for biological_run in 0..=run {
                            let true_start = if prefix { 70 + run - biological_run } else { 0 };
                            let true_end = if prefix {
                                150
                            } else {
                                80 - run + biological_run
                            };
                            assert!(core_start >= true_start && core_end <= true_end);
                        }
                        let value = serde_json::to_value(diagnostic.unwrap()).unwrap();
                        assert_eq!(
                            value["boundary_core_proposals"][0]["uncertainty_bases"],
                            run
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn boundary_core_merged_provenance_is_order_independent_and_ablatable() {
        let mut genome = crate::fragment::tests::body(8000);
        genome[2073..2080].fill(b'A');
        genome[2072] = b'C';
        let seq = [&genome[2000..2080], &vec![b'A'; 70]].concat();
        let idx = index(&genome);
        let original = candidate(
            2000,
            Strand::Plus,
            vec![CigarOp::Match(80), CigarOp::SoftClip(70)],
        );
        let alternative = candidate(
            2000,
            Strand::Plus,
            vec![CigarOp::Match(79), CigarOp::SoftClip(71)],
        );
        for order in [
            vec![original.clone(), alternative.clone()],
            vec![alternative.clone(), original.clone()],
        ] {
            let cached = cache(order);
            let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
            let result = aligner
                .recover_short_body_observed(&seq, &mate(2300, Strand::Minus), Some(&cached), None)
                .unwrap();
            assert_eq!(
                result.cigar,
                vec![CigarOp::Match(73), CigarOp::SoftClip(77)]
            );
            assert_eq!(result.boundary_core_uncertainty, Some((73, 7)));
            let view = crate::bam::aln_view(&result);
            assert!(view.tags.iter().any(|tag| tag.0 == *b"XB"
                && matches!(&tag.1, esperanto_bamio::TagValue::Str(s)
                    if s=="query_start=73;length=7;assignment=unknown")));
            aligner.config.boundary_core = false;
            assert!(aligner
                .recover_short_body_observed(&seq, &mate(2300, Strand::Minus), Some(&cached), None)
                .is_none());
        }
        // An already qualifying full exact-body decision must take precedence
        // over any new partial proposal from another source path.
        let (genome, seq, p) = fixture();
        let idx = index(&genome);
        let cached = cache(vec![p]);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let mut diagnostic = None;
        assert!(aligner
            .recover_short_body_observed(
                &seq,
                &mate(2300, Strand::Minus),
                Some(&cached),
                Some(&mut diagnostic)
            )
            .is_some());
        assert!(diagnostic.unwrap().boundary_core_proposals.is_empty());
    }

    #[test]
    fn boundary_core_rejects_invalid_sources_and_unsupported_mates() {
        let mut genome = crate::fragment::tests::body(8000);
        genome[2076..2080].fill(b'A');
        genome[2075] = b'C';
        let seq = [&genome[2000..2080], &vec![b'A'; 70]].concat();
        let idx = index(&genome);
        let p = candidate(
            2000,
            Strand::Plus,
            vec![CigarOp::Match(80), CigarOp::SoftClip(70)],
        );
        for location in [10usize, 100] {
            let mut bad = seq.clone();
            bad[location] = b'N';
            assert!(boundary_core(&idx, &bad, &p).is_none());
        }
        let mut mismatch = seq.clone();
        mismatch[10] = if mismatch[10] == b'A' { b'C' } else { b'A' };
        assert!(boundary_core(&idx, &mismatch, &p).is_none());
        for cigar in [
            vec![
                CigarOp::Match(40),
                CigarOp::Ins(1),
                CigarOp::Match(39),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::Match(40),
                CigarOp::Del(1),
                CigarOp::Match(40),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::Match(40),
                CigarOp::RefSkip(20),
                CigarOp::Match(40),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::SoftClip(1),
                CigarOp::Match(79),
                CigarOp::SoftClip(70),
            ],
        ] {
            assert!(boundary_core(&idx, &seq, &candidate(2000, Strand::Plus, cigar)).is_none());
        }
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let cached = cache(vec![p]);
        for mut m in [
            mate(2300, Strand::Plus),
            mate(4000, Strand::Minus),
            mate(2300, Strand::Minus),
        ] {
            if m.pos == 2300 && m.strand == Strand::Minus {
                m.rescued = true;
            }
            assert!(aligner
                .recover_short_body_observed(&seq, &m, Some(&cached), None)
                .is_none());
        }
        for quality in [0, 19, 255] {
            let mut m = mate(2300, Strand::Minus);
            m.mapq_override = Some(quality);
            assert!(aligner
                .recover_short_body_observed(&seq, &m, Some(&cached), None)
                .is_none());
        }
        let run_index = index(&vec![b'A'; 8000]);
        assert!(boundary_core(
            &run_index,
            &vec![b'A'; 150],
            &candidate(
                2000,
                Strand::Plus,
                vec![CigarOp::Match(80), CigarOp::SoftClip(70)]
            )
        )
        .is_none());
    }

    #[test]
    fn boundary_core_nearby_omitted_repeat_refuses() {
        let mut genome = crate::fragment::tests::body(8000);
        genome[2076..2080].fill(b'A');
        genome[2075] = b'C';
        let body = genome[2000..2080].to_vec();
        genome[2500..2580].copy_from_slice(&body);
        let seq = [body.as_slice(), &vec![b'A'; 70]].concat();
        let idx = index(&genome);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        assert!(aligner
            .recover_short_body_observed(
                &seq,
                &mate(2300, Strand::Minus),
                Some(&cache(vec![candidate(
                    2000,
                    Strand::Plus,
                    vec![CigarOp::Match(80), CigarOp::SoftClip(70)]
                )])),
                None
            )
            .is_none());
    }

    #[test]
    fn boundary_core_distinct_supported_projections_refuse_in_either_order() {
        let mut genome = crate::fragment::tests::body(8000);
        let mut middle = genome[3000..3080].to_vec();
        middle[..4].fill(b'A');
        middle[4] = b'G';
        middle[75] = b'G';
        middle[76..].fill(b'C');
        let seq = [&vec![b'A'; 90][..], &middle, &vec![b'C'; 90]].concat();
        genome[2000..2260].copy_from_slice(&seq);
        genome[3000..3080].fill(b'G');
        let idx = index(&genome);
        let left = candidate(
            2000,
            Strand::Plus,
            vec![CigarOp::Match(170), CigarOp::SoftClip(90)],
        );
        let right = candidate(
            2090,
            Strand::Plus,
            vec![CigarOp::SoftClip(90), CigarOp::Match(170)],
        );
        // Each projection has independent indexed support, but they emit
        // different sets of columns. A shared locus cannot erase that choice.
        for p in [&left, &right] {
            let (mut core, _, _, _) = boundary_core(&idx, &seq, p).unwrap();
            let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
            aligner.rescore_strict(&seq, &mut core);
            assert!(indexed_support(&idx, &seq, &core).is_some());
        }
        for order in [
            vec![left.clone(), right.clone()],
            vec![right.clone(), left.clone()],
        ] {
            let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
            let mut diagnostic = None;
            assert!(aligner
                .recover_short_body_observed(
                    &seq,
                    &mate(2400, Strand::Minus),
                    Some(&cache(order)),
                    Some(&mut diagnostic),
                )
                .is_none());
            let pool = diagnostic.unwrap();
            assert_eq!(pool.boundary_core_proposals.len(), 2);
            let decision = pool.short_body.unwrap();
            assert_eq!(decision.decision, "ambiguous_paths");
            assert_eq!(decision.supported_candidates.len(), 2);
        }
    }

    #[test]
    fn boundary_core_effective_fragment_lifts_uncertainty_and_seed_coordinates() {
        let mut genome = crate::fragment::tests::body(8000);
        genome[2076..2080].fill(b'A');
        genome[2075] = b'C';
        let seq = [&genome[2000..2080], &vec![b'A'; 70]].concat();
        let raw = revcomp(&seq);
        let idx = index(&genome);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let mut diagnostic = None;
        let result = aligner
            .recover_short_body_observed(
                &raw,
                &mate(2300, Strand::Plus),
                Some(&cache(vec![candidate(
                    2000,
                    Strand::Minus,
                    vec![CigarOp::Match(80), CigarOp::SoftClip(70)],
                )])),
                Some(&mut diagnostic),
            )
            .unwrap();
        assert_eq!(result.boundary_core_uncertainty, Some((76, 4)));
        let record = crate::fastq::FastqRecord {
            name: b"core-fragment/1".to_vec(),
            seq: [raw.as_slice(), &vec![b'C'; 50]].concat(),
            qual: vec![b'I'; 200],
        };
        let mut out = Vec::new();
        crate::diagnostics::write_row_with_fragment(
            &mut out,
            &record,
            Some(1),
            true,
            crate::diagnostics::CandidatePool::capture(&idx, 150, &[], "independent"),
            diagnostic,
            Some(crate::fragment::FragmentEvidence {
                fragment_length: 150,
                overlap_bases: 100,
                viable_offsets: 1,
                adapter_family: "synthetic",
                adapter_bases_r1: 50,
                adapter_bases_r2: 50,
            }),
        )
        .unwrap();
        let row: serde_json::Value = serde_json::from_slice(&out).unwrap();
        assert_eq!(
            row["mate_attempt"]["boundary_core_proposals"][0]["uncertainty_query_start"],
            126
        );
        assert_eq!(
            row["mate_attempt"]["boundary_core_proposals"][0]["uncertainty_bases"],
            4
        );
    }

    #[test]
    fn short_body_tail_bases_prefix_suffix_both_strands_and_wrong_path() {
        let (mut genome, _, _) = fixture();
        for tail in b"ACGT" {
            // All tail bases and both ends receive a true nonmatching boundary.
            let different = *b"ACGT".iter().find(|b| *b != tail).unwrap();
            genome[2000] = different;
            genome[2079] = different;
            let idx = index(&genome);
            for prefix in [false, true] {
                let aligned = &genome[2000..2080];
                let seq = if prefix {
                    [&vec![*tail; 70], aligned].concat()
                } else {
                    [aligned, &vec![*tail; 70]].concat()
                };
                for strand in [Strand::Plus, Strand::Minus] {
                    let ops = if prefix {
                        vec![CigarOp::SoftClip(70), CigarOp::Match(80)]
                    } else {
                        vec![CigarOp::Match(80), CigarOp::SoftClip(70)]
                    };
                    let p = candidate(2000, strand, ops);
                    let support = indexed_support(&idx, &seq, &p).unwrap();
                    assert!(support.len() >= 2);
                    assert!(support.windows(2).all(|s| s[0].0 + 15 <= s[1].0));
                    let raw = if strand == Strand::Minus {
                        revcomp(&seq)
                    } else {
                        seq.clone()
                    };
                    let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
                    let result = aligner
                        .recover_short_body_observed(
                            &raw,
                            &mate(
                                2300,
                                if strand == Strand::Minus {
                                    Strand::Plus
                                } else {
                                    Strand::Minus
                                },
                            ),
                            Some(&cache(vec![p.clone()])),
                            None,
                        )
                        .unwrap();
                    assert_eq!((result.pos, result.strand), (2000, strand));
                    assert_eq!(crate::mapq::mapq_of(&result), 0);
                    assert!(result.rescued);
                    let mut wrong = p;
                    wrong.pos += 200;
                    assert!(indexed_support(&idx, &seq, &wrong).is_none());
                }
            }
        }
    }

    #[test]
    fn short_body_rejects_unobservable_boundaries_and_unsupported_geometry() {
        let (genome, seq, p) = fixture();
        let idx = index(&genome);
        let mut same = seq.clone();
        same[80..].fill(seq[79]);
        assert!(
            indexed_support(&idx, &same, &p).is_none(),
            "matching biological and incidental technical endpoints are unknown"
        );
        for (at, base) in [(80, b'N'), (100, b'N'), (100, seq[79]), (20, b'N')] {
            let mut query = seq.clone();
            query[at] = base;
            assert!(indexed_support(&idx, &query, &p).is_none());
        }
        for cigar in [
            vec![
                CigarOp::Match(40),
                CigarOp::Ins(1),
                CigarOp::Match(39),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::Match(40),
                CigarOp::Del(1),
                CigarOp::Match(40),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::Match(40),
                CigarOp::SoftClip(70),
                CigarOp::Match(40),
            ],
            vec![
                CigarOp::RefSkip(1),
                CigarOp::Match(80),
                CigarOp::SoftClip(70),
            ],
            vec![
                CigarOp::Match(80),
                CigarOp::SoftClip(70),
                CigarOp::RefSkip(1),
            ],
            vec![CigarOp::Match(0), CigarOp::Match(80), CigarOp::SoftClip(70)],
        ] {
            assert!(indexed_support(&idx, &seq, &candidate(2000, Strand::Plus, cigar)).is_none());
        }
        let short = [&genome[2000..2020], &vec![seq[80]; 130]].concat();
        assert!(indexed_support(
            &idx,
            &short,
            &candidate(
                2000,
                Strand::Plus,
                vec![CigarOp::Match(20), CigarOp::SoftClip(130)]
            )
        )
        .is_none());
        let primary = candidate(2000, Strand::Plus, vec![CigarOp::Match(80)]);
        assert!(indexed_support(&idx, &genome[2000..2080], &primary).is_none());
    }

    #[test]
    fn short_body_mate_confidence_span_boundary_and_independent_switch() {
        let (genome, seq, p) = fixture();
        let idx = index(&genome);
        let c = cache(vec![p]);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let good = mate(2850, Strand::Minus); // end3000 - start2000 =1000 inclusive
        assert!(aligner
            .recover_short_body_observed(&seq, &good, Some(&c), None)
            .is_some());
        let mut bad = Vec::new();
        bad.push(mate(2851, Strand::Minus));
        bad.push(mate(2300, Strand::Plus));
        let mut wrong_contig = good.clone();
        wrong_contig.contig = 1;
        bad.push(wrong_contig);
        for q in [0, 19, 255] {
            let mut m = good.clone();
            m.mapq_override = Some(q);
            bad.push(m);
        }
        let mut rescued = good.clone();
        rescued.rescued = true;
        bad.push(rescued);
        let mut missing_geometry = good.clone();
        missing_geometry.cigar.clear();
        bad.push(missing_geometry);
        for m in bad {
            assert!(aligner
                .recover_short_body_observed(&seq, &m, Some(&c), None)
                .is_none());
        }
        aligner.config.mate_short_body = false;
        assert!(aligner
            .recover_short_body_observed(&seq, &good, Some(&c), None)
            .is_none());
        assert!(
            aligner
                .align_read_with_mate(&seq, 0, 2300, Some(&c))
                .is_none(),
            "legacy coordinate-only mate API cannot grant this qualification"
        );
        assert!(!AlignConfig::default().mate_short_body);
        assert!(!AlignConfig::rna_default().exact_path_competition);
    }

    #[test]
    fn short_body_duplicate_index_and_identical_missing_source_limit() {
        let (mut genome, seq, p) = fixture();
        let first = index(&genome);
        assert!(indexed_support(&first, &seq, &p).is_some());
        // Same represented observations cannot identify an unrepresented source.
        assert_eq!(
            indexed_support(&first, &seq, &p),
            indexed_support(&first, &seq, &p)
        );
        let body = genome[2000..2080].to_vec();
        genome[5000..5080].copy_from_slice(&body);
        let duplicated = index(&genome);
        assert!(indexed_support(&duplicated, &seq, &p).is_none());
        let mut aligner = Aligner::new(&duplicated, AlignConfig::rna_default());
        let conditional = aligner
            .recover_short_body_observed(
                &seq,
                &mate(2300, Strand::Minus),
                Some(&cache(vec![p])),
                None,
            )
            .expect("remote represented duplicate is compatible with conditional support");
        assert_eq!(crate::mapq::mapq_of(&conditional), 0);
        // Identical observations also arise when the source is missing or the
        // mate is wrong. Conditional evidence cannot identify those worlds.
    }

    #[test]
    fn short_body_conditional_remote_repeat_support_both_strands_and_clip_ends() {
        let (mut genome, _, _) = fixture();
        // Same complex body occurs remotely; no globally unique body seed can
        // authorize recovery. The supplied mate distinguishes the two copies.
        let body = genome[2000..2080].to_vec();
        genome[5000..5080].copy_from_slice(&body);
        for prefix in [false, true] {
            let adjacent = if prefix { body[0] } else { body[79] };
            let tail = *b"ACGT".iter().find(|&&b| b != adjacent).unwrap();
            let oriented = if prefix {
                [&vec![tail; 70], body.as_slice()].concat()
            } else {
                [body.as_slice(), &vec![tail; 70]].concat()
            };
            let idx = index(&genome);
            for strand in [Strand::Plus, Strand::Minus] {
                let cigar = if prefix {
                    vec![CigarOp::SoftClip(70), CigarOp::Match(80)]
                } else {
                    vec![CigarOp::Match(80), CigarOp::SoftClip(70)]
                };
                let p = candidate(2000, strand, cigar);
                assert!(indexed_support(&idx, &oriented, &p).is_none());
                let seq = if strand == Strand::Minus {
                    revcomp(&oriented)
                } else {
                    oriented.clone()
                };
                let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
                let mut diagnostic = None;
                let result = aligner.recover_short_body_observed(
                    &seq,
                    &mate(2300, strand.flip()),
                    Some(&cache(vec![p])),
                    Some(&mut diagnostic),
                );
                assert!(
                    result.is_some(),
                    "remote repeat must be distinguishable by existing complete mate geometry"
                );
                let result = result.unwrap();
                assert_eq!(result.pos, 2000);
                assert_eq!(crate::mapq::mapq_of(&result), 0);
                assert!(result.rescued && result.short_body_qualified);
                let evidence = serde_json::to_value(diagnostic.unwrap()).unwrap();
                assert_eq!(
                    evidence["short_body"]["supported_candidates"][0]["support_scope"],
                    "mate_conditioned_index"
                );
            }
        }
    }

    #[test]
    fn short_body_conditional_nearby_omitted_repeat_still_refuses() {
        let (mut genome, seq, p) = fixture();
        let body = genome[2000..2080].to_vec();
        genome[2500..2580].copy_from_slice(&body);
        let idx = index(&genome);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        // The competing full path is absent from the cached pool. Its indexed
        // seed occurrences must still defeat conditional support.
        assert!(aligner
            .recover_short_body_observed(
                &seq,
                &mate(2300, Strand::Minus),
                Some(&cache(vec![p])),
                None,
            )
            .is_none());
    }

    #[test]
    fn short_body_conditional_whole_list_budget_and_over_cap_refusal() {
        let (base, seq, p) = fixture();
        let segment = [&base[2000..2080], &base[3000..3040]].concat();
        for copies in [4500, 5001] {
            let mut genome = base.clone();
            for _ in 0..copies {
                genome.extend_from_slice(&segment);
            }
            let idx = index(&genome);
            assert!(indexed_support(&idx, &seq, &p).is_none());
            let mut attempt = ShortBodySketchAttempt::default();
            let supported =
                mate_conditioned_support(&idx, &seq, &p, &mate(2300, Strand::Minus), &mut attempt);
            assert!(attempt.scanned_postings <= MID_ROUND_MAX_ANCHORS);
            assert!(attempt.omitted_seed_lists > 0);
            if copies == 5001 {
                assert!(supported.is_none());
                assert_eq!(attempt.scanned_postings, 0);
                assert!(attempt.over_cap_seed_lists > 0);
            } else {
                assert!(attempt.scanned_postings > 0);
                assert_eq!(attempt.scanned_postings % (copies + 1), 0);
                let mut again = ShortBodySketchAttempt::default();
                assert_eq!(
                    supported,
                    mate_conditioned_support(
                        &idx,
                        &seq,
                        &p,
                        &mate(2300, Strand::Minus),
                        &mut again,
                    ),
                );
                assert_eq!(
                    serde_json::to_value(attempt).unwrap(),
                    serde_json::to_value(again).unwrap(),
                );
            }
        }
    }

    #[test]
    fn short_body_conditional_diagnostics_lift_seeds_without_global_unique_claim() {
        let (mut genome, seq, mut p) = fixture();
        let body = genome[2000..2080].to_vec();
        genome[5000..5080].copy_from_slice(&body);
        p.strand = Strand::Minus;
        let idx = index(&genome);
        let raw = revcomp(&seq);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let mut diagnostic = None;
        assert!(aligner
            .recover_short_body_observed(
                &raw,
                &mate(2300, Strand::Plus),
                Some(&cache(vec![p])),
                Some(&mut diagnostic),
            )
            .is_some());
        let pool = diagnostic.unwrap();
        let supported = &pool.short_body.as_ref().unwrap().supported_candidates[0];
        assert!(supported.disjoint_unique_indexed_seeds.is_empty());
        let original = supported.disjoint_mate_conditioned_indexed_seeds.clone();
        assert!(original.len() >= 2);
        let read = crate::fastq::FastqRecord {
            name: b"conditional-lift/1".to_vec(),
            seq: [raw.as_slice(), &vec![b'C'; 50]].concat(),
            qual: vec![b'I'; 200],
        };
        let mut output = Vec::new();
        crate::diagnostics::write_row_with_fragment(
            &mut output,
            &read,
            Some(1),
            true,
            crate::diagnostics::CandidatePool::capture(&idx, 150, &[], "independent"),
            Some(pool),
            Some(crate::fragment::FragmentEvidence {
                fragment_length: 150,
                overlap_bases: 100,
                viable_offsets: 1,
                adapter_family: "synthetic",
                adapter_bases_r1: 50,
                adapter_bases_r2: 50,
            }),
        )
        .unwrap();
        let row: serde_json::Value = serde_json::from_slice(&output).unwrap();
        let supported = &row["mate_attempt"]["short_body"]["supported_candidates"][0];
        assert_eq!(supported["support_scope"], "mate_conditioned_index");
        assert_eq!(
            supported["disjoint_unique_indexed_seeds"],
            serde_json::json!([])
        );
        let lifted: Vec<(u32, u64)> =
            serde_json::from_value(supported["disjoint_mate_conditioned_indexed_seeds"].clone())
                .unwrap();
        assert_eq!(
            lifted,
            original
                .into_iter()
                .map(|(q, key)| (q + 50, key))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn short_body_diagnostics_lift_minus_seed_offsets_with_effective_fragment() {
        let (genome, seq, mut p) = fixture();
        let idx = index(&genome);
        p.strand = Strand::Minus;
        let raw_effective = revcomp(&seq);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let mut diagnostic = None;
        assert!(aligner
            .recover_short_body_observed(
                &raw_effective,
                &mate(2300, Strand::Plus),
                Some(&cache(vec![p])),
                Some(&mut diagnostic)
            )
            .is_some());
        let pool = diagnostic.unwrap();
        let original = pool.short_body.as_ref().unwrap().supported_candidates[0]
            .disjoint_unique_indexed_seeds
            .clone();
        let record = crate::fastq::FastqRecord {
            name: b"effective/1".to_vec(),
            seq: [raw_effective.as_slice(), &vec![b'C'; 50]].concat(),
            qual: vec![b'I'; 200],
        };
        let independent = crate::diagnostics::CandidatePool::capture(&idx, 150, &[], "independent");
        let mut output = Vec::new();
        crate::diagnostics::write_row_with_fragment(
            &mut output,
            &record,
            Some(1),
            true,
            independent,
            Some(pool),
            Some(crate::fragment::FragmentEvidence {
                fragment_length: 150,
                overlap_bases: 100,
                viable_offsets: 1,
                adapter_family: "synthetic",
                adapter_bases_r1: 50,
                adapter_bases_r2: 50,
            }),
        )
        .unwrap();
        let row: serde_json::Value = serde_json::from_slice(&output).unwrap();
        assert_eq!(row["mate_attempt"]["candidates"][0]["cigar"], "50S80M70S");
        let lifted: Vec<(u32, u64)> = serde_json::from_value(
            row["mate_attempt"]["short_body"]["supported_candidates"][0]
                ["disjoint_unique_indexed_seeds"]
                .clone(),
        )
        .unwrap();
        assert_eq!(
            lifted,
            original
                .into_iter()
                .map(|(q, k)| (q + 50, k))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            row["mate_attempt"]["short_body"]["query_coordinate_frame"],
            "raw_alignment_query"
        );
    }

    #[test]
    fn short_body_spliced_competition_keeps_exact_paths_and_is_order_independent() {
        let (mut genome, _, _) = fixture();
        // Shared unique 50bp block supplies disjoint seeds. The duplicated 30bp
        // suffix leaves two distinct splice paths observationally unresolved.
        let suffix = genome[2200..2230].to_vec();
        genome[2400..2430].copy_from_slice(&suffix);
        let tail = *b"ACGT".iter().find(|&&b| b != suffix[29]).unwrap();
        let seq = [&genome[2000..2050], &suffix, &vec![tail; 70]].concat();
        let idx = index(&genome);
        let a = candidate(
            2000,
            Strand::Plus,
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(150),
                CigarOp::Match(30),
                CigarOp::SoftClip(70),
            ],
        );
        let b = candidate(
            2000,
            Strand::Plus,
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(350),
                CigarOp::Match(30),
                CigarOp::SoftClip(70),
            ],
        );
        assert!(indexed_support(&idx, &seq, &a).is_some());
        assert!(indexed_support(&idx, &seq, &b).is_some());
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        let m = mate(2500, Strand::Minus);
        for candidates in [vec![a.clone(), b.clone()], vec![b.clone(), a.clone()]] {
            let mut diag = None;
            assert!(aligner
                .recover_short_body_observed(&seq, &m, Some(&cache(candidates)), Some(&mut diag))
                .is_none());
            let decision = diag.unwrap().short_body.unwrap();
            assert_eq!(decision.decision, "ambiguous_paths");
            assert_eq!(decision.qualified_paths, 2);
        }
        let mut diag = None;
        assert!(aligner
            .recover_short_body_observed(
                &seq,
                &m,
                Some(&cache(vec![a.clone(), a])),
                Some(&mut diag)
            )
            .is_some());
        assert_eq!(diag.unwrap().short_body.unwrap().qualified_paths, 1);
    }
}
