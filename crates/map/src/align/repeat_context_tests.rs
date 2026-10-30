//! Native-index witnesses. All sequence is generated, never human-derived.
use super::*;

#[test]
fn larger_context_budgets_do_not_claim_unique_repeat_placement() {
    let (index, query) = fixture(true, true, 6000);
    for seq in [query.clone(), revcomp(&query)] {
        for (postings, candidates) in [(20_000, 8), (80_000, 8), (20_000, 32), (80_000, 32)] {
            let mut config = AlignConfig::rna_default();
            config.context_posting_budget = postings;
            config.context_proposal_budget = candidates;
            let mut aligner = Aligner::new(&index, config);
            let (result, cache) = aligner.align_read_full(&seq);
            let result = result.expect("repeat candidates must remain placeable");
            assert_eq!(crate::mapq::mapq_of(&result), 0);
            assert!(result.seed_search_limited);
            let budget = cache.unwrap().seed_budget.unwrap().repeat_context.unwrap();
            assert!(budget.scanned_postings <= postings);
            assert!(budget.retained_candidates <= candidates);
        }
    }
}

fn random(state: &mut u64, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            b"ACGT"[(*state & 3) as usize]
        })
        .collect()
}

pub(crate) fn fixture(editing: bool, duplicate: bool, copies: usize) -> (Index, Vec<u8>) {
    fixture_extra(editing, duplicate, copies, None)
}

fn fixture_extra(
    editing: bool,
    duplicate: bool,
    copies: usize,
    extra: Option<&[u8]>,
) -> (Index, Vec<u8>) {
    let mut state = 719;
    let mut body = random(&mut state, 112);
    for p in [92, 100, 107] {
        body[p] = b'A';
    }
    let mut query = body.clone();
    for p in [92, 100, 107] {
        query[p] = if editing { b'G' } else { b'C' };
    }
    let repeated = [vec![b'N'; 32], body[..100].to_vec(), vec![b'N'; 44]]
        .concat()
        .repeat(copies);
    let background = random(&mut state, 1_000_000);
    let target = [
        random(&mut state, 120),
        body.clone(),
        random(&mut state, 120),
    ]
    .concat();
    let mut fasta = format!(
        ">repeats\n{}\n>background\n{}\n>target\n{}\n",
        std::str::from_utf8(&repeated).unwrap(),
        std::str::from_utf8(&background).unwrap(),
        std::str::from_utf8(&target).unwrap()
    );
    if duplicate {
        let second = [random(&mut state, 120), body, random(&mut state, 120)].concat();
        fasta.push_str(&format!(
            ">indistinguishable\n{}\n",
            std::str::from_utf8(&second).unwrap()
        ));
    }
    if let Some(extra) = extra {
        fasta.push_str(&format!(
            ">strong\n{}\n",
            std::str::from_utf8(extra).unwrap()
        ));
    }
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(fasta.as_bytes()).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    assert!(
        index.freq_cutoff < 5000,
        "native frequency cutoff, never forced"
    );
    (index, query)
}

fn recovery(editing: bool, minus: bool) {
    let (index, seq) = fixture(editing, false, 6000);
    let seq = if minus { revcomp(&seq) } else { seq };
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let a = aligner
        .align_read(&seq)
        .expect("ordinary repetitive candidate must survive occurrence filtering");
    assert_eq!(
        (a.contig, a.pos, a.strand, a.score),
        (
            2,
            120,
            if minus { Strand::Minus } else { Strand::Plus },
            206
        )
    );
    assert_eq!(a.cigar, vec![CigarOp::Match(112)]);
}

#[test]
fn ordinary_repeat_context_ag_plus() {
    recovery(true, false);
}
#[test]
fn ordinary_repeat_context_ag_minus() {
    recovery(true, true);
}
#[test]
fn ordinary_repeat_context_transversion_plus() {
    recovery(false, false);
}
#[test]
fn ordinary_repeat_context_transversion_minus() {
    recovery(false, true);
}

#[test]
fn ordinary_repeat_context_reports_truthful_budget_and_cannot_raise_mapq() {
    let (index, seq) = fixture(false, false, 6000);
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let (a, cache) = aligner.align_read_full(&seq);
    let mut a = a.unwrap();
    let cache = cache.unwrap();
    let b = cache.seed_budget.unwrap();
    let context = b.repeat_context.unwrap();
    assert!(b.limited() && cache.seed_search_limited() && a.seed_search_limited);
    assert_eq!(
        a.origin.label(),
        "repeat_context",
        "additional contexts must not pretend to be chain/DP-derived candidates"
    );
    assert_eq!(a.chain_score, 0);
    assert_eq!(b.eligible_anchors, 0);
    assert_eq!(b.retained_anchors, 0);
    assert_eq!(b.omitted_seeds, 0);
    assert_eq!(context.query_length, 112);
    assert_eq!(
        context.eligible_seed_lists,
        context.scanned_seed_lists + context.omitted_seed_lists
    );
    assert!(context.scanned_postings <= 20000 && context.scanned_postings > 5000);
    assert!(context.unique_intervals <= context.scanned_postings);
    assert!(context.retained_candidates <= 8);
    assert_eq!(
        context.positive_candidates,
        context.retained_candidates + context.omitted_candidates
    );
    a.mapq_override = Some(60);
    a.rescued = true;
    assert_eq!(a.mapq_of_checked(), 0);
    let pool = cache.diagnostic_pool(&index, 112, "test");
    let json = serde_json::to_value(pool).unwrap();
    assert_eq!(json["seed_budget"]["repeat_context"]["query_length"], 112);
    assert!(json["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .any(|p| p["reference_name"] == "target" && p["common_score"] == 206));
}

#[test]
fn ordinary_repeat_context_equal_copies_remain_ambiguous() {
    let (index, seq) = fixture(false, true, 6000);
    let a = Aligner::new(&index, AlignConfig::rna_default())
        .align_read(&seq)
        .unwrap();
    assert_eq!(a.score, 206);
    assert_eq!(a.second_score, 206);
    assert!(a.n_cands >= 2);
    assert_eq!(a.mapq_of_checked(), 0);
}

#[test]
fn ordinary_repeat_context_oversize_empty_cache_survives_mate_fallback() {
    let (index, seq) = fixture(false, false, 21000);
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let (a, cache) = aligner.align_read_full(&seq);
    assert!(a.is_none());
    let cache = cache.unwrap();
    assert!(cache.placements.is_empty());
    assert!(cache.seed_search_limited());
    let b = cache.seed_budget.unwrap().repeat_context.unwrap();
    assert_eq!(b.scanned_postings, 0);
    assert_eq!(b.scanned_seed_lists, 0);
    assert!(b.omitted_seed_lists > 0);
    let mut observed = None;
    let a = aligner
        .align_read_with_mate_observed(&seq, 2, 120, Some(&cache), Some(&mut observed))
        .expect("mate window may still propose the stronger locus");
    assert_eq!((a.contig, a.pos), (2, 120));
    assert!(a.seed_search_limited);
    assert_eq!(a.mapq_of_checked(), 0);
    assert!(observed
        .unwrap()
        .seed_budget
        .unwrap()
        .repeat_context
        .is_some());
}

#[test]
fn ordinary_repeat_context_ablation_and_other_modes_preserve_parent_search() {
    let (index, seq) = fixture(false, false, 6000);
    for (rna, dir, enabled) in [
        (true, EditDir::Off, false),
        (false, EditDir::Off, true),
        (true, EditDir::Ag, true),
        (true, EditDir::Tc, true),
    ] {
        let mut cfg = AlignConfig::rna_default();
        cfg.rna = rna;
        cfg.extend.edit_dir = dir;
        cfg.repeat_context = enabled;
        let (a, cache) = Aligner::new(&index, cfg).align_read_full(&seq);
        if !enabled {
            assert!(a.is_none());
        }
        assert!(cache.is_none_or(|c| c.seed_budget.is_none_or(|b| b.repeat_context.is_none())));
    }
    let (index, body) = super::prefix_tests::fixture();
    let (_, cache) = Aligner::new(&index, AlignConfig::rna_default()).align_read_full(&body);
    assert!(
        cache
            .unwrap()
            .seed_budget
            .is_none_or(|b| b.repeat_context.is_none()),
        "perfect edge hit bypasses extra search"
    );
    let raw = [body, b"AGATCGGA".to_vec()].concat();
    let (_, cache) = Aligner::new(&index, AlignConfig::rna_default()).align_read_full(&raw);
    assert!(
        cache
            .unwrap()
            .seed_budget
            .is_none_or(|b| b.repeat_context.is_none()),
        "successful first chain bypasses extra search"
    );
}

#[test]
fn ordinary_repeat_context_reuse_matches_every_fresh_cached_candidate() {
    let (index, seq) = fixture(false, false, 6000);
    let reads = [
        seq.clone(),
        revcomp(&seq),
        vec![b'N'; 112],
        Vec::new(),
        index.reference.contigs[1].slice_ascii(200, 350),
        seq,
    ];
    let cfg = AlignConfig::rna_default();
    let mut reused = Aligner::new(&index, cfg);
    for q in reads.iter().chain(reads.iter().rev()) {
        let observed = reused.align_read_full(q);
        let expected = Aligner::new(&index, cfg).align_read_full(q);
        assert_eq!(format!("{:?}", observed.0), format!("{:?}", expected.0));
        let summarize = |c: Option<CachedPlacements>| {
            c.map(|c| {
                (
                    c.n_seeds,
                    c.genomic_any,
                    format!("{:?}", c.seed_budget),
                    c.placements
                        .iter()
                        .map(|p| {
                            format!(
                                "{:?}",
                                (
                                    p.origin,
                                    p.contig,
                                    p.pos,
                                    p.strand,
                                    &p.cigar,
                                    &p.junctions,
                                    p.chain_score,
                                    p.ext_score,
                                    &p.strict,
                                    p.ea_dir
                                )
                            )
                        })
                        .collect::<Vec<_>>(),
                )
            })
        };
        assert_eq!(summarize(observed.1), summarize(expected.1));
    }
}

#[test]
fn ordinary_repeat_context_candidates_do_not_veto_stronger_indel_or_splice() {
    let (_, seq) = fixture(false, false, 6000);
    for spliced in [false, true] {
        let extra = if spliced {
            [
                &seq[..56],
                b"GT".as_slice(),
                &vec![b'N'; 96],
                b"AG".as_slice(),
                &seq[56..],
            ]
            .concat()
        } else {
            [&seq[..56], &seq[57..]].concat()
        };
        let (index, q) = fixture_extra(false, false, 6000, Some(&extra));
        assert_eq!(q, seq);
        let (new, budget) =
            crate::repeat_context::generate(&index, &q, ExtendParams::default()).unwrap();
        assert!(new.iter().any(|p| p.contig == 2 && p.score == 206));
        let mut placements: Vec<_> = new
            .into_iter()
            .map(|p| {
                Placement::from_genomic(p.contig, p.pos, p.strand, p.cigar, Vec::new(), 0, p.score)
            })
            .collect();
        let (cigar, junctions) = if spliced {
            let junction = crate::gtf::Junction {
                contig: 3,
                start: 56,
                end: 156,
                minus_strand: false,
            };
            (
                vec![
                    CigarOp::Match(56),
                    CigarOp::RefSkip(100),
                    CigarOp::Match(56),
                ],
                vec![RefinedJunction {
                    junction,
                    signal: SpliceSignal::GtAg,
                    known_support: 1,
                }],
            )
        } else {
            (
                vec![CigarOp::Match(56), CigarOp::Ins(1), CigarOp::Match(55)],
                Vec::new(),
            )
        };
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        if spliced {
            aligner.jlib = Some(std::sync::Arc::new(JunctionLib::build([
                junctions[0].junction
            ])));
        }
        placements.push(Placement::from_genomic(
            3,
            0,
            Strand::Plus,
            cigar.clone(),
            junctions,
            70,
            0,
        ));
        for p in &mut placements {
            aligner.rescore_strict(&q, p);
        }
        dedup_by_footprint(&mut placements, q.len());
        sort_placements(&mut placements);
        let native = placements.iter().find(|p| p.contig == 3).unwrap();
        assert!(native.strict.score > 206);
        let a = aligner
            .arbitrate_with_budget(
                &q,
                placements,
                30,
                true,
                Some(crate::index::AnchorBudget {
                    repeat_context: Some(budget),
                    ..Default::default()
                }),
            )
            .unwrap();
        assert_eq!((a.contig, a.pos), (3, 0));
        assert_eq!(a.cigar, cigar);
        assert_eq!(a.mapq_of_checked(), 0);
    }
}

#[test]
fn ordinary_repeat_context_hyper_rescoring_does_not_add_genomic_dp_reextension() {
    let genome = vec![b'A'; 112];
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(
            &[
                b">context\n".as_slice(),
                &genome,
                b"\n>trigger\n",
                &genome,
                b"\n",
            ]
            .concat(),
        )
        .unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    let query = [vec![b'A'; 90], vec![b'G'; 22]].concat();
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let cigar = vec![CigarOp::Match(90), CigarOp::SoftClip(22)];
    let mut context =
        Placement::from_genomic(0, 0, Strand::Plus, cigar.clone(), Vec::new(), 0, 180);
    context.origin = crate::mapq::Origin::RepeatContext;
    let mut trigger = Placement::from_genomic(
        1,
        0,
        Strand::Plus,
        vec![CigarOp::Match(112)],
        Vec::new(),
        0,
        92,
    );
    trigger.origin = crate::mapq::Origin::L1;
    let mut placements = vec![context, trigger];
    for p in &mut placements {
        aligner.rescore_strict(&query, p);
    }
    let rc = revcomp(&query);
    aligner.hyper_rearbitrate(&query, &rc, &mut placements);
    let context = placements.iter().find(|p| p.contig == 0).unwrap();
    assert_eq!(context.cigar, cigar);
    assert_eq!(context.ea_dir, EditDir::Ag);
    assert_eq!(context.strict.score, 180);
    assert_eq!(context.chain_score, 0);
}
