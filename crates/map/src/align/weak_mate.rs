//! Default-off experiment; a weak mate proposes a locus but never prunes it.
use super::*;

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> (
        Index<'static>,
        Vec<u8>,
        Placement,
        ReadAlignment,
        ReadAlignment,
    ) {
        let mut genome = crate::fragment::tests::body(8000);
        let seq = genome[2000..2150].to_vec();
        genome[5000..5150].copy_from_slice(&seq);
        for i in [20, 50, 80, 110] {
            genome[5000 + i] = if seq[i] == b'A' { b'C' } else { b'A' };
        }
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(
                &[b">c\n".as_slice(), &genome, b"\n>remote\n", &seq, b"\n"].concat(),
            )
            .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let original = Placement::from_genomic(
            0,
            5000,
            Strand::Plus,
            vec![CigarOp::Match(150)],
            Vec::new(),
            30,
            276,
        );
        let old = ReadAlignment {
            contig: 0,
            pos: 5000,
            strand: Strand::Plus,
            cigar: original.cigar.clone(),
            score: 276,
            mapq_override: Some(1),
            ..Default::default()
        };
        let mate = ReadAlignment {
            contig: 0,
            pos: 2050,
            strand: Strand::Minus,
            cigar: vec![CigarOp::Match(150)],
            score: 300,
            mapq_override: Some(0),
            seed_search_limited: true,
            ..Default::default()
        };
        (index, seq, original, old, mate)
    }

    #[test]
    fn weak_exact_mate_proposes_without_claiming_confidence() {
        let (idx, seq, p, old, mate) = setup();
        let mut mate_body = Vec::new();
        idx.reference.contigs[0].decode_into(2050, 2200, &mut mate_body);
        let mate_seq = revcomp(&mate_body);
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        aligner.config.weak_mate_proposals = true;
        let mut p = p;
        aligner.rescore_strict(&seq, &mut p);
        let cache = CachedPlacements {
            placements: vec![p],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let before = cache.placements[0].cigar.clone();
        let mut diagnostic = None;
        let observed = aligner.propose_from_weak_mate(
            &seq,
            &old,
            &mate_seq,
            &mate,
            &cache,
            Some(&mut diagnostic),
        );
        assert!(
            observed.is_some(),
            "an exact weak mate must be usable as a proposal, not an anchor"
        );
        let result = observed.unwrap();
        assert_eq!(
            (result.pos, result.strand, result.score),
            (2000, Strand::Plus, 300)
        );
        assert_eq!(crate::mapq::mapq_of(&result), 0);
        assert!(result.seed_search_limited && result.rescued);
        let pool = diagnostic.unwrap();
        assert_eq!(pool.stage, "weak_mate_proposal");
        assert!(pool.candidates.iter().any(|p| p.reference_start0 == 5000));
        assert_eq!(cache.placements.len(), 1);
        assert_eq!(cache.placements[0].cigar, before);
        assert_eq!(cache.placements[0].pos, 5000);
    }

    fn mate_sequence(index: &Index) -> Vec<u8> {
        let mut body = Vec::new();
        index.reference.contigs[0].decode_into(2050, 2200, &mut body);
        revcomp(&body)
    }

    #[test]
    fn remote_tie_is_not_pruned_by_pair_geometry() {
        let (idx, seq, p, old, mate) = setup();
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        aligner.config.weak_mate_proposals = true;
        let remote = Placement::from_genomic(
            1,
            0,
            Strand::Plus,
            vec![CigarOp::Match(150)],
            Vec::new(),
            1,
            300,
        );
        let cache = CachedPlacements {
            placements: vec![p, remote],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let mut diagnostic = None;
        assert!(aligner
            .propose_from_weak_mate(
                &seq,
                &old,
                &mate_sequence(&idx),
                &mate,
                &cache,
                Some(&mut diagnostic)
            )
            .is_none());
        let pool = diagnostic.unwrap();
        assert!(pool.candidates.iter().any(|p| p.reference_name == "remote"));
        assert!(pool.candidates.iter().any(|p| p.reference_start0 == 2000));
        assert_eq!(cache.placements.len(), 2);
    }

    #[test]
    fn stronger_remote_path_cannot_be_deleted() {
        let (idx, mut seq, p, mut old, mate) = setup();
        seq[120..].fill(b'T');
        let mut genome = Vec::new();
        idx.reference.contigs[0].decode_into(0, 8000, &mut genome);
        genome[2120..2150].fill(b'C');
        let reference = crate::fasta::parse_fasta_bytes(
            &[b">c\n".as_slice(), &genome, b"\n>remote\n", &seq, b"\n"].concat(),
        )
        .unwrap();
        let idx = Index::build(reference, SeedParams { k: 15, w: 5 });
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        aligner.config.weak_mate_proposals = true;
        old.score = 1;
        let remote = Placement::from_genomic(
            1,
            0,
            Strand::Plus,
            vec![CigarOp::Match(150)],
            Vec::new(),
            1,
            300,
        );
        let cache = CachedPlacements {
            placements: vec![p, remote],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let mut diagnostic = None;
        assert!(aligner
            .propose_from_weak_mate(
                &seq,
                &old,
                &mate_sequence(&idx),
                &mate,
                &cache,
                Some(&mut diagnostic)
            )
            .is_none());
        let pool = diagnostic.unwrap();
        let remote = pool
            .candidates
            .iter()
            .find(|p| p.reference_name == "remote")
            .unwrap();
        let proposal = pool
            .candidates
            .iter()
            .find(|p| p.reference_start0 == 2000)
            .unwrap();
        assert!(remote.common_score > proposal.common_score);
    }

    #[test]
    fn default_off_and_invalid_proposers_refuse() {
        let (idx, seq, p, old, mate) = setup();
        let cache = CachedPlacements {
            placements: vec![p],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let original_mate_seq = mate_sequence(&idx);
        for case in 0..8 {
            let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
            aligner.config.weak_mate_proposals = case != 0;
            let mut mate = mate.clone();
            let mut old = old.clone();
            let mut mate_seq = original_mate_seq.clone();
            match case {
                0 => {}
                1 => {
                    mate.seed_search_limited = false;
                    mate.mapq_override = Some(20);
                }
                2 => mate.seed_search_limited = false,
                3 => mate.rescued = true,
                4 => old.mapq_override = Some(20),
                5 => mate_seq[50] = b'N',
                6 => mate.cigar = vec![CigarOp::Match(75), CigarOp::Ins(1), CigarOp::Match(74)],
                7 => aligner.config.rna = false,
                _ => unreachable!(),
            }
            assert!(
                aligner
                    .propose_from_weak_mate(&seq, &old, &mate_seq, &mate, &cache, None)
                    .is_none(),
                "case={case}"
            );
        }
    }

    #[test]
    fn wrong_orientation_and_no_strict_improvement_refuse() {
        let (idx, seq, p, old, mate) = setup();
        let cache = CachedPlacements {
            placements: vec![p],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        for case in 0..2 {
            let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
            aligner.config.weak_mate_proposals = true;
            let mut old = old.clone();
            let mut mate = mate.clone();
            let mut mseq = mate_sequence(&idx);
            if case == 0 {
                mate.strand = Strand::Plus;
                mseq = revcomp(&mseq);
            } else {
                old.score = 300;
            }
            assert!(aligner
                .propose_from_weak_mate(&seq, &old, &mseq, &mate, &cache, None)
                .is_none());
        }
    }

    #[test]
    fn reverse_strand_proposal_uses_reference_oriented_exact_geometry() {
        let (idx, seq, mut p, mut old, mut mate) = setup();
        p.strand = Strand::Minus;
        old.strand = Strand::Minus;
        mate.strand = Strand::Plus;
        let cache = CachedPlacements {
            placements: vec![p],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        aligner.config.weak_mate_proposals = true;
        let result = aligner
            .propose_from_weak_mate(
                &revcomp(&seq),
                &old,
                &revcomp(&mate_sequence(&idx)),
                &mate,
                &cache,
                None,
            )
            .unwrap();
        assert_eq!(
            (result.pos, result.strand, result.score),
            (2000, Strand::Minus, 300)
        );
        assert_eq!(crate::mapq::mapq_of(&result), 0);
    }

    #[test]
    fn existing_editing_domain_is_not_reinterpreted_as_an_ordinary_search() {
        let (idx, seq, mut p, old, mate) = setup();
        p.ea_dir = EditDir::Ag;
        let cache = CachedPlacements {
            placements: vec![p],
            n_seeds: 5,
            genomic_any: true,
            seed_budget: None,
        };
        let mut aligner = Aligner::new(&idx, AlignConfig::rna_default());
        aligner.config.weak_mate_proposals = true;
        assert!(aligner
            .propose_from_weak_mate(&seq, &old, &mate_sequence(&idx), &mate, &cache, None,)
            .is_none());
    }
}

impl Aligner<'_> {
    pub(crate) fn propose_from_weak_mate(
        &mut self,
        seq: &[u8],
        old: &ReadAlignment,
        mate_seq: &[u8],
        mate: &ReadAlignment,
        cache: &CachedPlacements,
        diagnostic: Option<&mut Option<crate::diagnostics::CandidatePool>>,
    ) -> Option<ReadAlignment> {
        if !self.config.rna
            || !self.config.weak_mate_proposals
            || old.rescued
            || mate.rescued
            || old.score <= 0
            || mate.score <= 0
            || crate::mapq::mapq_of(old) >= crate::mapq::PILEUP_GATE
            || crate::mapq::mapq_of(mate) != 0
            || !mate.seed_search_limited
            || cache.placements.iter().any(|p| p.ea_dir.is_on())
            || !exact_ordinary_path(self, mate_seq, mate)
        {
            return None;
        }
        let mut proposal = self.mate_window_placement(seq, mate.contig, mate.pos)?;
        let query = if proposal.strand == Strand::Minus {
            revcomp(seq)
        } else {
            seq.to_vec()
        };
        self.rescore_strict(&query, &mut proposal);
        let geometry = ReadAlignment {
            contig: proposal.contig,
            pos: proposal.pos,
            strand: proposal.strand,
            cigar: proposal.cigar.clone(),
            ..Default::default()
        };
        if !proposal.primary_eligible(seq.len())
            || !exact_ordinary_path(self, seq, &geometry)
            || crate::pair::relate(&geometry, mate) != crate::pair::PairRelation::Proper
        {
            return None;
        }
        let mut placements = cache.placements.clone();
        placements.push(Placement::from_genomic(
            old.contig,
            old.pos,
            old.strand,
            old.cigar.clone(),
            old.junctions.clone(),
            old.chain_score,
            old.score,
        ));
        placements.push(proposal.clone());
        for p in &mut placements {
            let query = if p.strand == Strand::Minus {
                revcomp(seq)
            } else {
                seq.to_vec()
            };
            self.rescore_strict(&query, p);
        }
        placements.sort_by_key(|p| std::cmp::Reverse(p.strict.score));
        if let Some(slot) = diagnostic {
            let mut pool = crate::diagnostics::CandidatePool::capture(
                self.index,
                seq.len(),
                &placements,
                "weak_mate_proposal",
            );
            pool.seed_budget = cache.seed_budget;
            *slot = Some(pool);
        }
        // A mate is an extra generator, not a license to delete a remote rival.
        if placements
            .iter()
            .any(|p| !same_path(p, &proposal) && p.strict.score >= proposal.strict.score)
            || proposal.strict.score <= old.score
        {
            return None;
        }
        let mut result =
            self.arbitrate_with_budget(seq, placements, cache.n_seeds, true, cache.seed_budget)?;
        if result.contig != proposal.contig
            || result.pos != proposal.pos
            || result.strand != proposal.strand
            || result.cigar != proposal.cigar
        {
            return None;
        }
        result.rescued = true;
        result.seed_search_limited = true;
        crate::mapq::force_mapq0(&mut result);
        Some(result)
    }
}

fn same_path(a: &Placement, b: &Placement) -> bool {
    a.contig == b.contig && a.pos == b.pos && a.strand == b.strand && a.cigar == b.cigar
}

fn exact_ordinary_path(aligner: &mut Aligner<'_>, sequence: &[u8], path: &ReadAlignment) -> bool {
    if sequence.is_empty()
        || sequence.len() > 256
        || path.cigar.is_empty()
        || aligner
            .index
            .reference
            .contigs
            .get(path.contig as usize)
            .is_none()
    {
        return false;
    }
    let query = if path.strand == Strand::Minus {
        revcomp(sequence)
    } else {
        sequence.to_vec()
    };
    let reference = &aligner.index.reference.contigs[path.contig as usize];
    let (mut q, mut rf) = (0usize, path.pos);
    let mut counts = [0usize; 4];
    for (i, op) in path.cigar.iter().enumerate() {
        match *op {
            CigarOp::Match(n) if n > 0 => {
                let Some(end) = rf.checked_add(n) else {
                    return false;
                };
                let Some(qe) = q.checked_add(n as usize) else {
                    return false;
                };
                if end > reference.len || qe > query.len() {
                    return false;
                }
                for &base in &query[q..qe] {
                    let Some(code) = b"ACGT".iter().position(|&b| b == base.to_ascii_uppercase())
                    else {
                        return false;
                    };
                    if reference.base(rf).to_ascii() != base.to_ascii_uppercase() {
                        return false;
                    }
                    counts[code] += 1;
                    rf += 1;
                }
                q = qe;
            }
            CigarOp::SoftClip(n) if n > 0 && (i == 0 || i + 1 == path.cigar.len()) => {
                let Some(end) = q.checked_add(n as usize) else {
                    return false;
                };
                if end > query.len() {
                    return false;
                }
                q = end;
            }
            _ => return false,
        }
    }
    let aligned: usize = counts.iter().sum();
    q == query.len()
        && aligned >= (query.len() as f64 * REFUSE_MIN_ALIGNED_FRAC).ceil() as usize
        && counts.iter().filter(|&&n| n > 0).count() >= 3
        && !counts.iter().any(|&n| n * 5 > aligned * 4)
}
