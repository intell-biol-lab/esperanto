//! A pair may establish a sparse substitution scoring domain, but never a
//! biological editing call or a calibrated location confidence.
use super::*;
use crate::diagnostics::PairSubstitutionDecision;

#[derive(Clone, Debug)]
pub(crate) struct PairWitness {
    dir: EditDir,
    positions: Vec<u32>,
}

fn complement(dir: EditDir) -> EditDir {
    match dir {
        EditDir::Ag => EditDir::Tc,
        EditDir::Tc => EditDir::Ag,
        EditDir::Off => EditDir::Off,
    }
}

/// An M/S path in reference-oriented query coordinates. Refuse internal clips,
/// gaps, empty operations, overflows and incomplete query consumption.
fn columns(
    index: &Index,
    contig: u32,
    pos: u32,
    cigar: &[CigarOp],
    len: usize,
) -> Option<Vec<Option<u32>>> {
    let reference = index.reference.contigs.get(contig as usize)?;
    let mut out = Vec::with_capacity(len);
    let mut rf = pos;
    for (i, op) in cigar.iter().enumerate() {
        match *op {
            CigarOp::Match(n) if n > 0 => {
                let end = rf.checked_add(n)?;
                if end > reference.len || out.len().checked_add(n as usize)? > len {
                    return None;
                }
                out.extend((rf..end).map(Some));
                rf = end;
            }
            CigarOp::SoftClip(n) if n > 0 && (i == 0 || i + 1 == cigar.len()) => {
                if out.len().checked_add(n as usize)? > len {
                    return None;
                }
                out.resize(out.len() + n as usize, None);
            }
            _ => return None,
        }
    }
    (out.len() == len && out.iter().any(Option::is_some)).then_some(out)
}

fn oriented(seq: &[u8], strand: Strand) -> Vec<u8> {
    if strand == Strand::Minus {
        revcomp(seq)
    } else {
        seq.iter().map(u8::to_ascii_uppercase).collect()
    }
}

fn canonical(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

fn direction(reference: u8, query: u8) -> Option<EditDir> {
    match (reference, query) {
        (b'A', b'G') => Some(EditDir::Ag),
        (b'T', b'C') => Some(EditDir::Tc),
        _ => None,
    }
}

/// All aligned mismatches must have one direction; only shared high quality
/// observations count as the witness. Input qualities are sequenced Phred+33.
fn observations(
    index: &Index,
    seq: &[u8],
    qual: &[u8],
    old: &ReadAlignment,
) -> Option<(EditDir, Vec<(u32, u8, u8)>)> {
    if seq.len() != qual.len() {
        return None;
    }
    let path = columns(index, old.contig, old.pos, &old.cigar, seq.len())?;
    let query = oriented(seq, old.strand);
    let reference = &index.reference.contigs[old.contig as usize];
    let mut dir = None;
    let mut changes = Vec::new();
    for (q, position) in path.iter().enumerate() {
        let Some(rf) = *position else { continue };
        let (a, b) = (reference.base(rf).to_ascii(), query[q]);
        if !canonical(a) || !canonical(b) {
            return None;
        }
        if a == b {
            continue;
        }
        let d = direction(a, b)?;
        if dir.is_some_and(|previous| previous != d) {
            return None;
        }
        dir = Some(d);
        let raw_q = if old.strand == Strand::Minus {
            seq.len() - 1 - q
        } else {
            q
        };
        changes.push((rf, b, qual[raw_q].checked_sub(33)?));
    }
    Some((dir?, changes))
}

pub(crate) fn witness(
    index: &Index,
    seqs: [&[u8]; 2],
    quals: [&[u8]; 2],
    old: [&ReadAlignment; 2],
) -> Option<PairWitness> {
    if old.iter().any(|a| {
        a.rescued
            || a.seed_search_limited
            || !(crate::mapq::PILEUP_GATE..255).contains(&crate::mapq::mapq_of(a))
    }) || old[0].contig != old[1].contig
        || old[0].strand == old[1].strand
        || crate::pair::relate(old[0], old[1]) != crate::pair::PairRelation::Proper
        || !old
            .iter()
            .any(|a| a.cigar.iter().any(|op| matches!(op, CigarOp::SoftClip(_))))
    {
        return None;
    }
    let (d1, a) = observations(index, seqs[0], quals[0], old[0])?;
    let (d2, b) = observations(index, seqs[1], quals[1], old[1])?;
    if d1 != d2 {
        return None;
    }
    let mut positions: Vec<_> = a
        .iter()
        .filter(|&&(rf, base, quality)| {
            quality >= 20 && b.iter().any(|&(r, v, q)| r == rf && v == base && q >= 20)
        })
        .map(|&(rf, _, _)| rf)
        .collect();
    positions.sort_unstable();
    positions.dedup();
    (positions.len() >= 2).then_some(PairWitness { dir: d1, positions })
}

fn preserves_and_extends(
    index: &Index,
    seq: &[u8],
    old: &ReadAlignment,
    candidate: &Placement,
    dir: EditDir,
) -> bool {
    if candidate.contig != old.contig || candidate.strand != old.strand {
        return false;
    }
    let Some(before) = columns(index, old.contig, old.pos, &old.cigar, seq.len()) else {
        return false;
    };
    let Some(after) = columns(
        index,
        candidate.contig,
        candidate.pos,
        &candidate.cigar,
        seq.len(),
    ) else {
        return false;
    };
    let query = oriented(seq, old.strand);
    let reference = &index.reference.contigs[old.contig as usize];
    let mut gain = false;
    for (q, (a, b)) in before.iter().zip(&after).enumerate() {
        if a.is_some() && a != b {
            return false;
        }
        if a.is_none() {
            if let Some(rf) = *b {
                let (r, v) = (reference.base(rf).to_ascii(), query[q]);
                if !canonical(r) || !canonical(v) || (r != v && direction(r, v) != Some(dir)) {
                    return false;
                }
                gain = true;
            }
        }
    }
    gain
}

impl Aligner<'_> {
    /// Candidate competition is independent for each mate and conditioned on
    /// one witness computed from the original pair snapshot. Never bootstrap
    /// the second mate from a first-mate extension.
    pub(crate) fn extend_pair_substitution_observed(
        &mut self,
        seq: &[u8],
        old: &ReadAlignment,
        cached: &CachedPlacements,
        evidence: &PairWitness,
        diagnostic: Option<&mut Option<crate::diagnostics::CandidatePool>>,
    ) -> Option<ReadAlignment> {
        if !self.config.rna || !self.config.pair_substitution || cached.seed_search_limited() {
            return None;
        }
        let mut original = Placement::from_genomic(
            old.contig,
            old.pos,
            old.strand,
            old.cigar.clone(),
            old.junctions.clone(),
            old.chain_score,
            old.score,
        );
        original.origin = old.origin;
        let old_query = oriented(seq, old.strand);
        self.rescore_with(&old_query, &mut original, evidence.dir);
        let original_score = original.strict.score;
        let mut placements = cached.placements.clone();
        placements.push(original);
        // Preserve original proposals alongside extended copies. Every retained
        // candidate uses this same read-oriented domain, including competitors.
        let extended: Vec<_> = placements
            .iter()
            .map(|p| {
                let mut p = p.clone();
                let dir = if p.strand == old.strand {
                    evidence.dir
                } else {
                    complement(evidence.dir)
                };
                self.reextend_with(&oriented(seq, p.strand), &mut p, dir);
                p
            })
            .collect();
        placements.extend(extended);
        for p in &mut placements {
            let dir = if p.strand == old.strand {
                evidence.dir
            } else {
                complement(evidence.dir)
            };
            self.rescore_with(&oriented(seq, p.strand), p, dir);
        }
        dedup_by_footprint_with_options(&mut placements, seq.len(), false);
        sort_placements(&mut placements);
        let accepted = placements
            .iter()
            .find(|p| p.primary_eligible(seq.len()))
            .is_some_and(|p| {
                p.strict.score > original_score
                    && preserves_and_extends(self.index, seq, old, p, evidence.dir)
            });
        if let Some(slot) = diagnostic {
            let mut pool = crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "pair_substitution",
            );
            pool.seed_budget = cached.seed_budget;
            pool.pair_substitution = Some(PairSubstitutionDecision {
                decision: if accepted {
                    "extended_mapq0"
                } else {
                    "no_safe_score_gain"
                },
                reference_direction: if evidence.dir == EditDir::Ag {
                    "ag"
                } else {
                    "tc"
                },
                shared_reference_positions0: evidence.positions.clone(),
                original_conditioned_score: original_score,
            });
            *slot = Some(pool);
        }
        if !accepted {
            return None;
        }
        let mut result = self.arbitrate_with_budget(
            seq,
            placements,
            cached.n_seeds,
            cached.genomic_any,
            cached.seed_budget,
        )?;
        result.rescued = true;
        result.pair_substitution_qualified = true;
        crate::mapq::force_mapq0(&mut result);
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(dir: EditDir, rival: Option<Strand>) -> (Index, [Vec<u8>; 2], [ReadAlignment; 2]) {
        let mut genome = crate::fragment::tests::body(8000);
        let mut query = genome[2000..2131].to_vec();
        let (reference, observed) = if dir == EditDir::Ag {
            (b'A', b'G')
        } else {
            (b'T', b'C')
        };
        for q in [1, 17, 50, 51] {
            genome[2000 + q] = reference;
            query[q] = observed;
        }
        if let Some(strand) = rival {
            genome[5000..5131].copy_from_slice(&oriented(&query, strand));
        }
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let old = [Strand::Plus, Strand::Minus].map(|strand| ReadAlignment {
            contig: 0,
            pos: 2002,
            strand,
            cigar: vec![CigarOp::SoftClip(2), CigarOp::Match(129)],
            score: 240,
            mapq_override: Some(60),
            ..Default::default()
        });
        (index, [query.clone(), revcomp(&query)], old)
    }

    fn get_witness(
        index: &Index,
        seqs: &[Vec<u8>; 2],
        old: &[ReadAlignment; 2],
    ) -> Option<PairWitness> {
        witness(
            index,
            [&seqs[0], &seqs[1]],
            [&vec![b'J'; 131], &vec![b'J'; 131]],
            [&old[0], &old[1]],
        )
    }

    fn cached(old: &ReadAlignment, rival: Option<Strand>) -> CachedPlacements {
        let mut placements = vec![Placement::from_genomic(
            old.contig,
            old.pos,
            old.strand,
            old.cigar.clone(),
            Vec::new(),
            30,
            240,
        )];
        if let Some(strand) = rival {
            placements.push(Placement::from_genomic(
                0,
                5000,
                strand,
                vec![CigarOp::Match(131)],
                Vec::new(),
                30,
                262,
            ));
        }
        CachedPlacements {
            placements,
            n_seeds: 10,
            genomic_any: true,
            seed_budget: None,
        }
    }

    #[test]
    fn witness_uses_both_raw_qualities_and_reference_orientations() {
        for dir in [EditDir::Ag, EditDir::Tc] {
            let (index, seqs, old) = fixture(dir, None);
            let w = get_witness(&index, &seqs, &old).unwrap();
            assert_eq!(w.dir, dir);
            assert_eq!(w.positions, vec![2017, 2050, 2051]);
            assert_eq!(
                get_witness(
                    &index,
                    &[seqs[1].clone(), seqs[0].clone()],
                    &[old[1].clone(), old[0].clone()]
                )
                .unwrap()
                .positions,
                w.positions,
            );
            let high = vec![b'J'; 131];
            let mut low = high.clone();
            for q in [17, 50, 51] {
                low[130 - q] = 33 + 19;
            }
            assert!(witness(
                &index,
                [&seqs[0], &seqs[1]],
                [&high, &low],
                [&old[0], &old[1]]
            )
            .is_none());
            low[130 - 17] = 33 + 20;
            assert!(witness(
                &index,
                [&seqs[0], &seqs[1]],
                [&high, &low],
                [&old[0], &old[1]]
            )
            .is_none());
            low[130 - 50] = 33 + 20;
            assert_eq!(
                witness(
                    &index,
                    [&seqs[0], &seqs[1]],
                    [&high, &low],
                    [&old[0], &old[1]]
                )
                .unwrap()
                .positions,
                vec![2017, 2050]
            );
            assert!(witness(
                &index,
                [&seqs[0], &seqs[1]],
                [&high[..130], &high],
                [&old[0], &old[1]]
            )
            .is_none());
        }
    }

    #[test]
    fn witness_refuses_weak_limited_rescued_or_incompatible_originals() {
        let (index, seqs, old) = fixture(EditDir::Tc, None);
        for mate in 0..2 {
            for case in 0..9 {
                let mut modified = old.clone();
                let other_strand = modified[1 - mate].strand;
                let m = &mut modified[mate];
                match case {
                    0 => m.mapq_override = Some(19),
                    1 => m.mapq_override = Some(255),
                    2 => m.rescued = true,
                    3 => m.seed_search_limited = true,
                    4 => m.contig = 1,
                    5 => m.strand = other_strand,
                    6 => m.pos = 6000,
                    7 => m.cigar = vec![CigarOp::Match(20), CigarOp::Ins(1), CigarOp::Match(110)],
                    _ => {
                        m.cigar = vec![
                            CigarOp::Match(20),
                            CigarOp::RefSkip(30),
                            CigarOp::Match(111),
                        ]
                    }
                }
                assert!(
                    get_witness(&index, &seqs, &modified).is_none(),
                    "mate={mate} case={case}"
                );
            }
        }
        let mut full = old.clone();
        for a in &mut full {
            a.pos = 2000;
            a.cigar = vec![CigarOp::Match(131)];
        }
        assert!(get_witness(&index, &seqs, &full).is_none());
    }

    #[test]
    fn witness_refuses_ordinary_mixed_unknown_and_unshared_changes() {
        let (index, seqs, old) = fixture(EditDir::Tc, None);
        for base in [b'A', b'G', b'N'] {
            let mut modified = seqs.clone();
            modified[0][17] = base; // reference T: ordinary or unknown
            assert!(get_witness(&index, &modified, &old).is_none());
        }
        let mut mixed = seqs.clone();
        // Manufacture a mismatching non-T reference column with the other
        // canonical editing direction when possible, otherwise an ordinary one.
        let rf = index.reference.contigs[0].base(2020).to_ascii();
        mixed[0][20] = if rf == b'A' { b'G' } else { b'A' };
        if mixed[0][20] == rf {
            mixed[0][20] = b'T';
        }
        assert!(get_witness(&index, &mixed, &old).is_none());
        let mut single = seqs.clone();
        for q in [50, 51] {
            single[1][130 - q] = b'A';
        } // restore mate to reference T
        assert!(get_witness(&index, &single, &old).is_none());
        let mut none = seqs.clone();
        for q in [17, 50, 51] {
            none[0][q] = b'T';
            none[1][130 - q] = b'A';
        }
        assert!(get_witness(&index, &none, &old).is_none());
    }

    #[test]
    fn transaction_extends_each_original_strand_with_truthful_score_and_mapq0() {
        for dir in [EditDir::Ag, EditDir::Tc] {
            let (index, seqs, old) = fixture(dir, None);
            let evidence = get_witness(&index, &seqs, &old).unwrap();
            let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
            for mate in 0..2 {
                let mut diagnostic = None;
                let result = aligner
                    .extend_pair_substitution_observed(
                        &seqs[mate],
                        &old[mate],
                        &cached(&old[mate], None),
                        &evidence,
                        Some(&mut diagnostic),
                    )
                    .unwrap();
                assert_eq!(result.pos, 2000);
                assert_eq!(result.cigar, vec![CigarOp::Match(131)]);
                assert_eq!(result.strand, old[mate].strand);
                assert_eq!(result.score, 254);
                assert_eq!(result.ea_count(), 4);
                assert_eq!(crate::mapq::mapq_of(&result), 0);
                assert!(result.rescued && result.pair_substitution_qualified);
                let d = diagnostic.unwrap().pair_substitution.unwrap();
                assert_eq!(d.original_conditioned_score, 252);
                assert_eq!(d.decision, "extended_mapq0");
            }
        }
    }

    #[test]
    fn all_competitors_use_one_domain_and_a_better_rival_refuses_extension() {
        for dir in [EditDir::Ag, EditDir::Tc] {
            for rival in [Strand::Plus, Strand::Minus] {
                let (index, seqs, old) = fixture(dir, Some(rival));
                let evidence = get_witness(&index, &seqs, &old).unwrap();
                let mut cache = cached(&old[0], Some(rival));
                for reversed in [false, true] {
                    if reversed {
                        cache.placements.reverse();
                    }
                    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
                    let mut diagnostic = None;
                    assert!(aligner
                        .extend_pair_substitution_observed(
                            &seqs[0],
                            &old[0],
                            &cache,
                            &evidence,
                            Some(&mut diagnostic),
                        )
                        .is_none());
                    let pool = diagnostic.unwrap();
                    assert_eq!(
                        pool.pair_substitution.unwrap().decision,
                        "no_safe_score_gain"
                    );
                    assert!(pool.candidates.iter().all(|c| c.editing_score_mode
                        == if (c.strand == "+") == (dir == EditDir::Ag) {
                            "ag"
                        } else {
                            "tc"
                        }));
                }
            }
        }
    }

    #[test]
    fn acceptance_requires_added_valid_columns_and_preserves_every_old_column() {
        let (index, seqs, old) = fixture(EditDir::Tc, None);
        let mut p = cached(&old[0], None).placements.remove(0);
        p.pos = 2000;
        p.cigar = vec![CigarOp::Match(131)];
        assert!(preserves_and_extends(
            &index,
            &seqs[0],
            &old[0],
            &p,
            EditDir::Tc
        ));
        for change in 0..5 {
            let mut bad = p.clone();
            match change {
                0 => bad.pos += 1,
                1 => bad.strand = Strand::Minus,
                2 => bad.cigar = vec![CigarOp::Match(60), CigarOp::Ins(1), CigarOp::Match(70)],
                3 => bad.cigar = old[0].cigar.clone(),
                _ => bad.cigar = vec![CigarOp::Match(60), CigarOp::SoftClip(1), CigarOp::Match(70)],
            }
            assert!(!preserves_and_extends(
                &index,
                &seqs[0],
                &old[0],
                &bad,
                EditDir::Tc
            ));
        }
        for base in [b'A', b'G', b'N'] {
            let mut bad = seqs[0].clone();
            bad[1] = base;
            assert!(!preserves_and_extends(
                &index,
                &bad,
                &old[0],
                &p,
                EditDir::Tc
            ));
        }
        assert!(columns(&index, 99, 0, &[CigarOp::Match(131)], 131).is_none());
        assert!(columns(&index, 0, u32::MAX, &[CigarOp::Match(131)], 131).is_none());
    }

    #[test]
    fn neutral_clipped_substitutions_cannot_gain_by_coverage_tie() {
        let (index, seqs, mut old) = fixture(EditDir::Tc, None);
        let evidence = get_witness(&index, &seqs, &old).unwrap();
        // Drop the preceding exact q0. The only clipped base is now the
        // known q1 T>C, whose addition scores zero in this domain.
        let query = seqs[0][1..].to_vec();
        for a in &mut old {
            a.cigar = vec![CigarOp::SoftClip(1), CigarOp::Match(129)];
        }
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        let mut diagnostic = None;
        assert!(aligner
            .extend_pair_substitution_observed(
                &query,
                &old[0],
                &cached(&old[0], None),
                &evidence,
                Some(&mut diagnostic),
            )
            .is_none());
        assert_eq!(
            diagnostic
                .unwrap()
                .pair_substitution
                .unwrap()
                .original_conditioned_score,
            252
        );
    }

    #[test]
    fn gap_rival_and_disabled_or_limited_transactions_preserve_original() {
        let (initial, seqs, old) = fixture(EditDir::Tc, None);
        let mut genome = initial.reference.contigs[0].slice_ascii(0, 8000);
        genome[5000..5065].copy_from_slice(&seqs[0][..65]);
        genome[5065] = b'A';
        genome[5066..5132].copy_from_slice(&seqs[0][65..]);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let evidence = get_witness(&index, &seqs, &old).unwrap();
        let mut cache = cached(&old[0], None);
        cache.placements.push(Placement::from_genomic(
            0,
            5000,
            Strand::Plus,
            vec![CigarOp::Match(65), CigarOp::Del(1), CigarOp::Match(66)],
            Vec::new(),
            30,
            256,
        ));
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        let mut diagnostic = None;
        assert!(aligner
            .extend_pair_substitution_observed(
                &seqs[0],
                &old[0],
                &cache,
                &evidence,
                Some(&mut diagnostic),
            )
            .is_none());
        assert_eq!(diagnostic.unwrap().candidates[0].common_score, 256);
        for dna in [false, true] {
            let mut cfg = AlignConfig::rna_default();
            if dna {
                cfg.rna = false;
            } else {
                cfg.pair_substitution = false;
            }
            let mut a = Aligner::new(&index, cfg);
            assert!(a
                .extend_pair_substitution_observed(
                    &seqs[0],
                    &old[0],
                    &cached(&old[0], None),
                    &evidence,
                    None,
                )
                .is_none());
        }
        let mut limited = cached(&old[0], None);
        limited.seed_budget = Some(crate::index::AnchorBudget {
            eligible_anchors: 100,
            retained_anchors: 50,
            omitted_seeds: 1,
            repeat_context: None,
        });
        assert!(aligner
            .extend_pair_substitution_observed(&seqs[0], &old[0], &limited, &evidence, None,)
            .is_none());
    }
}
