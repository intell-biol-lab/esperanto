//! Observable qualification of coverage-ineligible paired RNA candidates.
//! Indexed sketch support is not exhaustive uniqueness; output is always MAPQ0.
use super::*;
use crate::diagnostics::{ShortBodyCandidate, ShortBodyDecision};
use crate::index::unpack_pos;

/// Exact, complex body with one observably separated canonical terminal tail.
/// Query is already in candidate alignment orientation.
fn indexed_support(index: &Index, query: &[u8], p: &Placement) -> Option<Vec<(u32, u64)>> {
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
    let mut supported = Vec::new();
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
                    let hits = index.query(seed.kmer);
                    if hits.count != 1 {
                        continue;
                    }
                    let (contig, position, strand) = unpack_pos(hits.positions[0]);
                    if contig == p.contig && position == r + seed.pos && strand == seed.strand {
                        supported.push((q as u32 + seed.pos, seed.kmer));
                    }
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
    supported.sort_unstable();
    supported.dedup();
    let mut end = 0;
    let mut disjoint = Vec::new();
    for (start, key) in supported {
        if start >= end {
            disjoint.push((start, key));
            end = start + index.params.k;
        }
    }
    (disjoint.len() >= 2).then_some(disjoint)
}

fn same_path(a: &Placement, b: &Placement) -> bool {
    a.contig == b.contig && a.pos == b.pos && a.strand == b.strand && a.cigar == b.cigar
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
        let mut qualifying: Vec<(usize, Vec<(u32, u64)>)> = Vec::new();
        for (rank, p) in placements.iter().enumerate() {
            if p.contig != mate.contig || p.strand == mate.strand {
                continue;
            }
            let query = if p.strand == Strand::Minus { &rc } else { seq };
            let Some(seeds) = indexed_support(self.index, query, p) else {
                continue;
            };
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
            if !qualifying
                .iter()
                .any(|(i, _)| same_path(&placements[*i], p))
            {
                qualifying.push((rank, seeds));
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
                    .map(|(rank, seeds)| ShortBodyCandidate {
                        rank_in_pool: *rank,
                        disjoint_unique_indexed_seeds: seeds.clone(),
                    })
                    .collect(),
            });
            *slot = Some(pool);
        }
        self.scratch.rc = rc;
        if qualifying.len() != 1 {
            return None;
        }
        let (rank, seeds) = &qualifying[0];
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

    fn index(genome: &[u8]) -> Index {
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
        assert!(aligner
            .recover_short_body_observed(
                &seq,
                &mate(2300, Strand::Minus),
                Some(&cache(vec![p])),
                None
            )
            .is_none());
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
