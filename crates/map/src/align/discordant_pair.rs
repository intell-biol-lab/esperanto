//! Conditional use of already scored loci for discordant ordinary pairs.
use super::*;

fn full_match(cigar: &[CigarOp], length: usize) -> bool {
    !cigar.is_empty() && cigar == [CigarOp::Match(length as u32)]
}

/// Disjoint, inward-facing ordinary bodies provide separate reference support.
/// Overlapping bodies cannot supply the independent evidence this proposal needs.
fn compatible(a: &ReadAlignment, b: &ReadAlignment) -> bool {
    if a.contig != b.contig || a.strand == b.strand {
        return false;
    }
    let (plus, minus) = if a.strand == Strand::Plus {
        (a, b)
    } else {
        (b, a)
    };
    crate::pair::ref_end(plus) <= minus.pos
        && crate::pair::transcript_span(a, b) <= crate::pair::PROPER_SPAN_MAX
}

impl Aligner<'_> {
    pub(crate) fn propose_discordant_pair(
        &mut self,
        seq: &[u8],
        old: &ReadAlignment,
        mate_seq: &[u8],
        mate: &ReadAlignment,
        cache: &CachedPlacements,
        mate_cache: &CachedPlacements,
        diagnostic: Option<&mut Option<crate::diagnostics::CandidatePool>>,
    ) -> Option<ReadAlignment> {
        if !self.config.rna
            || !self.config.discordant_pair_arbitration
            || old.rescued
            || mate.rescued
            || old.seed_search_limited
            || mate.seed_search_limited
            || cache.seed_search_limited()
            || mate_cache.seed_search_limited()
            || !full_match(&old.cigar, seq.len())
            || !full_match(&mate.cigar, mate_seq.len())
            || crate::mapq::mapq_of(old) < crate::mapq::PILEUP_GATE
            || crate::mapq::mapq_of(old) == 255
            || crate::mapq::mapq_of(mate) < crate::mapq::PILEUP_GATE
            || crate::mapq::mapq_of(mate) == 255
            || crate::pair::relate(old, mate) == crate::pair::PairRelation::Proper
            || !weak_mate::exact_ordinary_path(self, mate_seq, mate)
            || cache.placements.iter().any(|p| p.ea_dir.is_on())
            || mate_cache.placements.iter().any(|p| p.ea_dir.is_on())
        {
            return None;
        }
        // Only an unchanged independently scored target is eligible. Earlier
        // rescues/extensions must not borrow the original candidate confidence.
        let original = cache.placements.iter().find(|p| {
            p.contig == old.contig
                && p.pos == old.pos
                && p.strand == old.strand
                && p.cigar == old.cigar
                && p.strict.score == old.score
        })?;
        let placements: Vec<_> = cache
            .placements
            .iter()
            .filter(|p| {
                let path = ReadAlignment {
                    contig: p.contig,
                    pos: p.pos,
                    strand: p.strand,
                    cigar: p.cigar.clone(),
                    ..Default::default()
                };
                compatible(&path, mate)
            })
            .cloned()
            .collect();
        let best = placements.iter().find(|p| p.primary_eligible(seq.len()))?;
        if !full_match(&best.cigar, seq.len())
            || best.origin != crate::mapq::Origin::Genomic
            || is_decoy_contig(self.index.reference, best.contig)
            || best.strict.score < original.strict.score - crate::mapq::UNIQUE_MARGIN
        {
            return None;
        }
        // A mate at the wrong homolog can drag an exact target. Compare both
        // fragment hypotheses in the unchanged ordinary score domain, including
        // eligible clipped reciprocal paths. A competing pair must trail by
        // at least the existing UNIQUE_MARGIN; weak remote paths alone cannot
        // veto a separated pair. This never promotes the reciprocal path.
        let proposed_pair_score = best.strict.score + mate.score;
        if mate_cache.placements.iter().any(|p| {
            p.primary_eligible(mate_seq.len())
                && old.score + p.strict.score > proposed_pair_score - crate::mapq::UNIQUE_MARGIN
                && compatible(
                    &ReadAlignment {
                        contig: p.contig,
                        pos: p.pos,
                        strand: p.strand,
                        cigar: p.cigar.clone(),
                        ..Default::default()
                    },
                    old,
                )
        }) {
            return None;
        }
        // Preserve within-fragment competitors in the ordinary scorer. A local
        // tie is not unique and cannot acquire pileup confidence by this route.
        let mut result = self.arbitrate_with_budget(
            seq,
            placements.clone(),
            cache.n_seeds,
            cache.genomic_any,
            cache.seed_budget,
        )?;
        if crate::mapq::mapq_of(&result) < crate::mapq::PILEUP_GATE {
            return None;
        }
        result.rescued = true; // existing ceiling is PILEUP_GATE, never above it
        if let Some(slot) = diagnostic {
            let mut pool = crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "discordant_pair_arbitration",
            );
            pool.seed_budget = cache.seed_budget;
            *slot = Some(pool);
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn discordant_pair_refuses_reciprocal_support_even_below_proposal_score_margin() {
        let mut genome = crate::fragment::tests::body(8000);
        let read = genome[1000..1100].to_vec();
        genome[5000..5100].copy_from_slice(&read);
        genome[5045] = if genome[5045] == b'A' { b'C' } else { b'A' };
        let mate_body = genome[1300..1400].to_vec();
        genome[5300..5400].copy_from_slice(&mate_body);
        for i in 5396..5400 {
            genome[i] = if genome[i] == b'A' { b'C' } else { b'A' };
        }
        let mate_seq = revcomp(&genome[5300..5400]);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let mut config = AlignConfig::rna_default();
        config.discordant_pair_arbitration = true;
        let mut aligner = Aligner::new(&index, config);
        let mut targets = vec![
            placement(0, 1000, &[CigarOp::Match(100)]),
            placement(0, 5000, &[CigarOp::Match(100)]),
        ];
        let mut mates = vec![
            placement(0, 5300, &[CigarOp::Match(100)]),
            placement(0, 1300, &[CigarOp::Match(96), CigarOp::SoftClip(4)]),
        ];
        for p in &mut targets {
            aligner.rescore_strict(&read, p);
        }
        for p in &mut mates {
            p.strand = Strand::Minus;
            aligner.rescore_strict(&revcomp(&mate_seq), p);
        }
        let old = aligner.arbitrate(&read, targets.clone(), 20, true).unwrap();
        let mate = aligner
            .arbitrate(&mate_seq, mates.clone(), 20, true)
            .unwrap();
        assert_eq!((old.score, old.second_score), (200, 194));
        assert_eq!((mate.score, mate.second_score), (200, 192));
        assert_eq!(crate::mapq::mapq_of(&mate), crate::mapq::PILEUP_GATE);
        let cache = CachedPlacements {
            placements: targets,
            n_seeds: 20,
            genomic_any: true,
            seed_budget: None,
        };
        let mate_cache = CachedPlacements {
            placements: mates,
            n_seeds: 20,
            genomic_any: true,
            seed_budget: None,
        };
        assert!(
            aligner
                .propose_discordant_pair(&read, &old, &mate_seq, &mate, &cache, &mate_cache, None)
                .is_none(),
            "a clipped reciprocal mate hypothesis must veto sacrificing an exact target, even though its score deficit exceeds the one-substitution proposal margin"
        );
    }

    #[test]
    fn discordant_pair_keeps_separated_joint_support_despite_weak_reciprocal_path() {
        let mut genome = crate::fragment::tests::body(8000);
        let mut read = genome[1000..1100].to_vec();
        read[45] = if read[45] == b'A' { b'G' } else { b'A' };
        genome[5000..5100].copy_from_slice(&read);
        let mate_body = genome[1300..1400].to_vec();
        genome[5300..5400].copy_from_slice(&mate_body);
        for b in &mut genome[5396..5400] {
            *b = if *b == b'A' { b'C' } else { b'A' };
        }
        let mate_seq = revcomp(&mate_body);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let mut config = AlignConfig::rna_default();
        config.discordant_pair_arbitration = true;
        let mut aligner = Aligner::new(&index, config);
        let mut targets = vec![
            placement(0, 5000, &[CigarOp::Match(100)]),
            placement(0, 1000, &[CigarOp::Match(100)]),
        ];
        let mut mates = vec![
            placement(0, 1300, &[CigarOp::Match(100)]),
            placement(0, 5300, &[CigarOp::Match(91), CigarOp::SoftClip(9)]),
        ];
        for p in &mut targets {
            aligner.rescore_strict(&read, p);
        }
        for p in &mut mates {
            p.strand = Strand::Minus;
            aligner.rescore_strict(&mate_body, p);
        }
        let old = aligner.arbitrate(&read, targets.clone(), 20, true).unwrap();
        let mate = aligner
            .arbitrate(&mate_seq, mates.clone(), 20, true)
            .unwrap();
        assert_eq!((old.score, old.second_score), (200, 194));
        assert_eq!((mate.score, mate.second_score), (200, 182));
        let cache = CachedPlacements {
            placements: targets,
            n_seeds: 20,
            genomic_any: true,
            seed_budget: None,
        };
        let mate_cache = CachedPlacements {
            placements: mates,
            n_seeds: 20,
            genomic_any: true,
            seed_budget: None,
        };
        let result = aligner
            .propose_discordant_pair(&read, &old, &mate_seq, &mate, &cache, &mate_cache, None)
            .expect("a weak reciprocal382 pair must not veto separated proposed394 support");
        assert_eq!(result.pos, 1000);
        assert_eq!(crate::mapq::mapq_of(&result), crate::mapq::PILEUP_GATE);
        for (aligned, accepted) in [(94, true), (96, false)] {
            let mut alternate = mate_cache.placements[1].clone();
            alternate.cigar = vec![CigarOp::Match(aligned), CigarOp::SoftClip(100 - aligned)];
            aligner.rescore_strict(&mate_body, &mut alternate);
            let mut alternatives = mate_cache.clone();
            alternatives.placements[1] = alternate;
            let supporting_mate = aligner
                .arbitrate(&mate_seq, alternatives.placements.clone(), 20, true)
                .unwrap();
            let result = aligner.propose_discordant_pair(
                &read,
                &old,
                &mate_seq,
                &supporting_mate,
                &cache,
                &alternatives,
                None,
            );
            assert_eq!(result.is_some(), accepted, "joint separation at {aligned}M");
        }
    }

    #[test]
    fn conditional_pair_requires_complete_exact_disjoint_mate_and_local_separation() {
        let mut genome = crate::fragment::tests::body(8000);
        let body = genome[1000..1100].to_vec();
        genome[1150..1250].copy_from_slice(&body);
        let mut read = body;
        read[45] = if read[45] == b'A' { b'G' } else { b'A' };
        genome[5000..5100].copy_from_slice(&read);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        for case in 0..17 {
            let mut config = AlignConfig::rna_default();
            config.discordant_pair_arbitration = case != 1;
            config.rna = case != 2;
            let mut aligner = Aligner::new(&index, config);
            let mut ps: Vec<_> = [5000, 1000]
                .into_iter()
                .map(|pos| {
                    Placement::from_genomic(
                        0,
                        pos,
                        Strand::Plus,
                        vec![CigarOp::Match(100)],
                        vec![],
                        200,
                        200,
                    )
                })
                .collect();
            for p in &mut ps {
                aligner.rescore_strict(&read, p);
            }
            let mut old = aligner.arbitrate(&read, ps.clone(), 20, true).unwrap();
            let mut mate_seq = revcomp(&genome[1300..1400]);
            let (mate, mate_cache) = aligner.align_read_full(&mate_seq);
            let mut mate = mate.unwrap();
            let mut mate_cache = mate_cache.unwrap();
            match case {
                3 => mate.rescued = true,
                4 => old.rescued = true,
                5 => old.seed_search_limited = true,
                6 => mate.seed_search_limited = true,
                7 => mate_seq[20] = if mate_seq[20] == b'A' { b'C' } else { b'A' },
                8 => {
                    mate.cigar = vec![CigarOp::Match(50), CigarOp::RefSkip(10), CigarOp::Match(50)]
                }
                9 => old.cigar = vec![CigarOp::Match(99), CigarOp::SoftClip(1)],
                10 => crate::mapq::force_mapq0(&mut mate),
                11 => {
                    let mut alternative = ps[1].clone();
                    alternative.pos = 1150;
                    ps.push(alternative);
                }
                12 => {
                    mate_seq = revcomp(&genome[1060..1160]);
                    mate = ReadAlignment {
                        contig: 0,
                        pos: 1060,
                        strand: Strand::Minus,
                        score: 200,
                        cigar: vec![CigarOp::Match(100)],
                        mapq_override: Some(60),
                        ..Default::default()
                    };
                }
                13 => ps[1].strict.score = ps[0].strict.score - crate::mapq::UNIQUE_MARGIN - 1,
                14 => crate::mapq::force_mapq0(&mut old),
                15 => mate_cache.placements[0].ea_dir = EditDir::Ag,
                16 => old.mapq_override = Some(255),
                _ => {}
            }
            let cache = CachedPlacements {
                placements: ps,
                n_seeds: 20,
                genomic_any: true,
                seed_budget: None,
            };
            let proposal = aligner.propose_discordant_pair(
                &read,
                &old,
                &mate_seq,
                &mate,
                &cache,
                &mate_cache,
                None,
            );
            assert_eq!(proposal.is_some(), case == 0, "guard case {case}");
        }
    }

    #[test]
    fn pair_support_geometry_rejects_overlap_outward_long_or_same_strand() {
        let a = ReadAlignment {
            contig: 0,
            pos: 1000,
            strand: Strand::Plus,
            cigar: vec![CigarOp::Match(100)],
            ..Default::default()
        };
        let b = ReadAlignment {
            contig: 0,
            pos: 1300,
            strand: Strand::Minus,
            cigar: vec![CigarOp::Match(100)],
            ..Default::default()
        };
        assert!(compatible(&a, &b));
        assert!(compatible(&b, &a));
        for case in 0..5 {
            let mut bad = b.clone();
            match case {
                0 => bad.pos = 1050,
                1 => bad.pos = 800,
                2 => bad.pos = 2000,
                3 => bad.strand = Strand::Plus,
                4 => bad.contig = 1,
                _ => unreachable!(),
            }
            assert!(!compatible(&a, &bad), "geometry case {case}");
        }
    }
}
