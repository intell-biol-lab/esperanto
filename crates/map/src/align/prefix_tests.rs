//! Synthetic read-boundary competition. No publication/human sequences.
use super::*;

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

pub(crate) fn fixture() -> (Index, Vec<u8>) {
    fixture_with_copies(0, 0)
}

fn fixture_with_copies(exact_copies: usize, interior_copies: usize) -> (Index, Vec<u8>) {
    let mut state = 1;
    let body = random(&mut state, 142);
    let mut decoy = body.clone();
    for p in [110, 140] {
        decoy[p] = b"ACGT"[(b"ACGT".iter().position(|&x| x == decoy[p]).unwrap() + 1) % 4];
    }
    let target = [
        random(&mut state, 120),
        body.clone(),
        random(&mut state, 120),
    ]
    .concat();
    let background = random(&mut state, 1_000_000);
    let repeats = [vec![b'N'; 32], decoy, vec![b'N'; 32]]
        .concat()
        .repeat(6000);
    let mut fasta = format!(
        ">decoys\n{}\n>background\n{}\n>target\n{}\n",
        std::str::from_utf8(&repeats).unwrap(),
        std::str::from_utf8(&background).unwrap(),
        std::str::from_utf8(&target).unwrap()
    );
    if exact_copies + interior_copies > 0 {
        let mut interior = body.clone();
        interior[140] =
            b"ACGT"[(b"ACGT".iter().position(|&x| x == interior[140]).unwrap() + 1) % 4];
        let copies = [
            [vec![b'N'; 32], body.clone(), vec![b'N'; 32]]
                .concat()
                .repeat(exact_copies),
            [vec![b'N'; 32], interior, vec![b'N'; 32]]
                .concat()
                .repeat(interior_copies),
        ]
        .concat();
        fasta.push_str(&format!(
            ">additional\n{}\n",
            std::str::from_utf8(&copies).unwrap()
        ));
    }
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(fasta.as_bytes()).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    // Native frequency policy, not a manually relaxed test cutoff.
    assert_eq!(index.freq_cutoff, 2);
    (index, body)
}

#[test]
fn query_edge_completion_recovers_exact_prefix_on_both_strands() {
    let (index, body) = fixture();
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let raw = [body.clone(), b"AGATCGGA".to_vec()].concat();
    let baseline = aligner
        .align_read(&raw)
        .expect("raw context supplies a body seed");
    assert_eq!(baseline.contig, 2);
    assert_eq!(baseline.pos, 120);
    assert_eq!(
        baseline.cigar,
        vec![CigarOp::Match(142), CigarOp::SoftClip(8)]
    );
    for (read, strand) in [
        (body.clone(), Strand::Plus),
        (revcomp(&body), Strand::Minus),
    ] {
        // The unchanged reference/read minimizer sketch has no chain. The
        // recovery must come from query-edge candidates, not a lower gate.
        let mins = crate::seed::minimizers(&read, aligner.config.seed);
        assert!(aligner
            .chain_round(&mins, 142, 15, MID_OCC_CAP)
            .is_none_or(|r| r.best().is_none()));
        let result = aligner
            .align_read(&read)
            .expect("an exact prefix must not need adapter context for its body seed");
        assert_eq!(result.contig, 2);
        assert_eq!(result.pos, 120);
        assert_eq!(result.strand, strand);
        assert_eq!(result.score, 284);
        assert_eq!(result.cigar, vec![CigarOp::Match(142)]);
        assert_eq!(result.mm_count, 0);
    }
}

#[test]
fn query_edge_completion_bypasses_successful_rounds_and_other_modes() {
    let (mut index, body) = fixture();
    let raw = [body.clone(), b"AGATCGGA".to_vec()].concat();
    let cfg = AlignConfig::rna_default();
    let original = crate::seed::minimizers(&raw, cfg.seed).len();
    let mut aligner = Aligner::new(&index, cfg);
    assert!(aligner
        .chain_round(
            &crate::seed::minimizers(&raw, cfg.seed),
            150,
            15,
            index.freq_cutoff
        )
        .unwrap()
        .best()
        .is_some());
    assert_eq!(aligner.genomic_candidates(&raw).n_seeds, original);

    // A successful bounded mid round also bypasses completion.
    index.freq_cutoff = 0;
    let mut aligner = Aligner::new(&index, cfg);
    assert!(aligner
        .chain_round(&crate::seed::minimizers(&raw, cfg.seed), 150, 15, 0)
        .is_none());
    assert!(aligner
        .bounded_chain_round(&crate::seed::minimizers(&raw, cfg.seed), 150, 15)
        .0
        .unwrap()
        .best()
        .is_some());
    assert_eq!(aligner.genomic_candidates(&raw).n_seeds, original);

    for (rna, edit_dir) in [
        (false, EditDir::Off),
        (true, EditDir::Ag),
        (true, EditDir::Tc),
    ] {
        let mut other = cfg;
        other.rna = rna;
        other.extend.edit_dir = edit_dir;
        let mut aligner = Aligner::new(&index, other);
        let original = crate::seed::minimizers(&body, other.seed).len();
        assert_eq!(aligner.genomic_candidates(&body).n_seeds, original);
    }
}

#[test]
fn query_edge_completion_does_not_force_negative_terminal_columns() {
    let (index, body) = fixture();
    let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
    let mut query = body.clone();
    // A terminal transversion has strictly negative terminal score,
    // rather than a neutral editing-type tie. The internal boundary seed
    // remains available, but is not authority to absorb these columns.
    query[0] = match body[0] {
        b'A' => b'C',
        b'C' => b'A',
        b'G' => b'T',
        _ => b'G',
    };
    for read in [query.clone(), revcomp(&query)] {
        let (result, cache) = aligner.align_read_full(&read);
        let original = crate::seed::minimizers(&read, aligner.config.seed).len();
        let cache = cache.expect("negative-tail witness must reach candidate arbitration");
        assert!(
            cache.n_seeds > original,
            "exercise the new edge fallback: {} <= {}",
            cache.n_seeds,
            original
        );
        let result = result.expect("a supported internal body remains alignable");
        assert_eq!(result.contig, 2);
        assert_eq!(result.pos, 121);
        assert_eq!(
            result.cigar,
            vec![CigarOp::SoftClip(1), CigarOp::Match(141)]
        );
        assert_eq!(result.score, 282);
        assert_eq!(result.mm_count, 0);
    }
}

#[test]
fn query_edge_completion_keeps_actual_overflow_search_limited() {
    let (index, body) = fixture_with_copies(999, 3999);
    for read in [body.clone(), revcomp(&body)] {
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default());
        let mins = crate::seed::minimizers(&read, aligner.config.seed);
        assert!(aligner
            .chain_round(&mins, 142, 15, index.freq_cutoff)
            .is_none_or(|r| r.best().is_none()));
        let (old, budget) = aligner.bounded_chain_round(&mins, 142, 15);
        assert!(
            old.is_none_or(|r| r.best().is_none()),
            "no old round may bypass fallback"
        );
        assert!(budget.limited());
        let (result, cache) = aligner.align_read_full(&read);
        let cache = cache.expect("limited candidate evidence must survive");
        assert!(cache.n_seeds > mins.len());
        let budget = cache.seed_budget.expect("actual overflow must be declared");
        assert!(budget.limited());
        assert!(budget.eligible_anchors > MID_ROUND_MAX_ANCHORS);
        assert!(budget.retained_anchors <= MID_ROUND_MAX_ANCHORS);
        let result = result.expect("budgeted edge seeds must retain a real repeat placement");
        assert_eq!(result.mapq_of_checked(), 0);
        assert_eq!(result.score, 284);
        assert_eq!(result.cigar, vec![CigarOp::Match(142)]);
        let reference = &index.reference.contigs[result.contig as usize];
        let mut observed = Vec::new();
        reference.decode_into(result.pos, result.pos + 142, &mut observed);
        assert_eq!(observed, body);
    }
}

#[test]
fn query_edge_completion_scratch_reuse_matches_fresh_aligners() {
    let (index, body) = fixture_with_copies(999, 3999);
    let reads = [
        [body.clone(), b"AGATCGGA".to_vec()].concat(),
        body.clone(),
        Vec::new(),
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_vec(),
        index.reference.contigs[1].slice_ascii(100, 250),
        revcomp(&body),
        b"ACGT".to_vec(),
    ];
    let cfg = AlignConfig::rna_default();
    let mut reused = Aligner::new(&index, cfg);
    for read in reads.iter().chain(reads.iter().rev()) {
        let observed = reused.align_read_full(read);
        let expected = Aligner::new(&index, cfg).align_read_full(read);
        assert_eq!(format!("{:?}", observed.0), format!("{:?}", expected.0));
        let summarize = |cache: Option<CachedPlacements>| {
            cache.map(|c| {
                let placements: Vec<_> = c
                    .placements
                    .iter()
                    .map(|p| {
                        let Placement {
                            origin,
                            contig,
                            pos,
                            strand,
                            cigar,
                            junctions,
                            chain_score,
                            ext_score,
                            strict,
                            ea_dir,
                        } = p;
                        format!(
                            "{:?}",
                            (
                                origin,
                                contig,
                                pos,
                                strand,
                                cigar,
                                junctions,
                                chain_score,
                                ext_score,
                                strict,
                                ea_dir
                            )
                        )
                    })
                    .collect();
                (
                    c.n_seeds,
                    c.genomic_any,
                    format!("{:?}", c.seed_budget),
                    placements,
                )
            })
        };
        assert_eq!(summarize(observed.1), summarize(expected.1));
    }
}
