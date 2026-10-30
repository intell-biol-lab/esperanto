//! Generated sequence-indistinguishable splice paths; no human truth inputs.
use super::*;
use crate::{fasta::parse_fasta_bytes, gtf::Junction, index::Index};

fn experimental_config() -> AlignConfig {
    let mut config = AlignConfig::rna_default();
    config.exact_path_competition = true;
    config
}

fn dedup_experimental(ps: &mut Vec<Placement>, read_len: usize) {
    dedup_by_footprint_with_options(ps, read_len, true);
}

fn random(n: usize) -> Vec<u8> {
    let mut state = 0x19eace11u64;
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn placement(pos: u32, cigar: Vec<CigarOp>, strand: Strand) -> Placement {
    Placement::from_genomic(0, pos, strand, cigar, Vec::new(), 150, 300)
}

fn terminal_fixture(
    n: usize,
    left: bool,
    distinct_base: bool,
) -> (Index<'static>, Vec<u8>, Arc<JunctionLib>, u32, Vec<CigarOp>) {
    let mut genome = random(2000);
    let query = genome[200..350].to_vec();
    let (pos, start, end, cigar) = if left {
        genome[100..100 + n].copy_from_slice(&query[..n]);
        (
            100,
            100 + n as u32,
            200 + n as u32,
            vec![
                CigarOp::Match(n as u32),
                CigarOp::RefSkip(100),
                CigarOp::Match(150 - n as u32),
            ],
        )
    } else {
        genome[1100..1100 + n].copy_from_slice(&query[150 - n..]);
        (
            200,
            350 - n as u32,
            1100,
            vec![
                CigarOp::Match(150 - n as u32),
                CigarOp::RefSkip(750 + n as u32),
                CigarOp::Match(n as u32),
            ],
        )
    };
    if distinct_base {
        let p = if left { 200 } else { 349 };
        genome[p] = if genome[p] == b'A' { b'C' } else { b'A' };
    }
    let fasta = [b">chrT\n".to_vec(), genome, b"\n".to_vec()].concat();
    let index = Index::build(
        parse_fasta_bytes(&fasta).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    let library = Arc::new(JunctionLib::build([Junction {
        contig: 0,
        start,
        end,
        minus_strand: false,
    }]));
    (index, query, library, pos, cigar)
}

#[test]
fn exact_paths_terminal_structures_survive_and_use_simple_ambiguous_primary() {
    for n in [1, 2, 20] {
        for left in [false, true] {
            for strand in [Strand::Plus, Strand::Minus] {
                let (index, query, library, pos, cigar) = terminal_fixture(n, left, false);
                let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library);
                let mut continuous = placement(200, vec![CigarOp::Match(150)], strand);
                let mut spliced = placement(pos, cigar, strand);
                aligner.rescore_strict(&query, &mut continuous);
                aligner.rescore_strict(&query, &mut spliced);
                assert_eq!(continuous.strict.score, 300);
                assert_eq!(spliced.strict.score, 306);
                let raw = if strand == Strand::Plus {
                    query.clone()
                } else {
                    revcomp(&query)
                };
                for mut paths in [
                    vec![spliced.clone(), continuous.clone(), continuous.clone()],
                    vec![continuous.clone(), spliced.clone()],
                ] {
                    dedup_experimental(&mut paths, raw.len());
                    assert_eq!(
                        paths.len(),
                        2,
                        "distinct exact mappings must survive footprint dedup"
                    );
                    let mut result = aligner.arbitrate(&raw, paths, 30, true).unwrap();
                    assert_eq!(result.cigar, vec![CigarOp::Match(150)]);
                    assert_eq!(crate::mapq::mapq_of(&result), 0);
                    result.mapq_override = Some(60);
                    assert_eq!(
                        crate::mapq::mapq_of(&result),
                        0,
                        "paired override cannot erase ambiguity"
                    );
                    assert_eq!(crate::mapq::margin(&result).2, "exact_path_tie");
                }
            }
        }
    }
}

#[test]
fn exact_paths_extra_junction_does_not_create_sequence_evidence() {
    let mut genome = random(2000);
    let query = [genome[200..275].to_vec(), genome[475..550].to_vec()].concat();
    genome[1100..1102].copy_from_slice(&query[148..]);
    let fasta = [b">chrT\n".to_vec(), genome, b"\n".to_vec()].concat();
    let index = Index::build(
        parse_fasta_bytes(&fasta).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    let library = Arc::new(JunctionLib::build([
        Junction {
            contig: 0,
            start: 275,
            end: 475,
            minus_strand: false,
        },
        Junction {
            contig: 0,
            start: 548,
            end: 1100,
            minus_strand: false,
        },
    ]));
    for strand in [Strand::Plus, Strand::Minus] {
        let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library.clone());
        let mut simple = placement(
            200,
            vec![
                CigarOp::Match(75),
                CigarOp::RefSkip(200),
                CigarOp::Match(75),
            ],
            strand,
        );
        let mut complex = placement(
            200,
            vec![
                CigarOp::Match(75),
                CigarOp::RefSkip(200),
                CigarOp::Match(73),
                CigarOp::RefSkip(552),
                CigarOp::Match(2),
            ],
            strand,
        );
        aligner.rescore_strict(&query, &mut simple);
        aligner.rescore_strict(&query, &mut complex);
        assert_eq!((simple.strict.score, complex.strict.score), (306, 312));
        let mut paths = vec![complex, simple.clone()];
        dedup_experimental(&mut paths, 150);
        assert_eq!(paths.len(), 2);
        let raw = if strand == Strand::Plus {
            query.clone()
        } else {
            revcomp(&query)
        };
        let result = aligner.arbitrate(&raw, paths, 30, true).unwrap();
        assert_eq!(result.cigar, simple.cigar);
        assert_eq!(crate::mapq::mapq_of(&result), 0);
    }
}

#[test]
fn exact_paths_sequence_supported_short_exon_still_wins() {
    for n in [1, 2] {
        for left in [false, true] {
            for strand in [Strand::Plus, Strand::Minus] {
                let (index, query, library, pos, cigar) = terminal_fixture(n, left, true);
                let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library);
                let mut continuous = placement(200, vec![CigarOp::Match(150)], strand);
                let mut spliced = placement(pos, cigar.clone(), strand);
                aligner.rescore_strict(&query, &mut continuous);
                aligner.rescore_strict(&query, &mut spliced);
                assert_eq!(continuous.strict.matches, 149);
                assert_eq!(spliced.strict.matches, 150);
                let mut paths = vec![continuous, spliced];
                dedup_experimental(&mut paths, 150);
                let raw = if strand == Strand::Plus {
                    query.clone()
                } else {
                    revcomp(&query)
                };
                let result = aligner.arbitrate(&raw, paths, 30, true).unwrap();
                assert_eq!(result.cigar, cigar);
                assert_eq!(crate::mapq::mapq_of(&result), 60);
                assert_ne!(crate::mapq::margin(&result).2, "exact_path_tie");
            }
        }
    }
}

#[test]
fn exact_paths_equivalent_cigar_splitting_is_not_ambiguity() {
    let (index, query, library, _, _) = terminal_fixture(1, false, false);
    let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library);
    let mut a = placement(200, vec![CigarOp::Match(150)], Strand::Plus);
    let mut b = placement(
        200,
        vec![CigarOp::Match(50), CigarOp::Match(100)],
        Strand::Plus,
    );
    aligner.rescore_strict(&query, &mut a);
    aligner.rescore_strict(&query, &mut b);
    let mut paths = vec![a, b];
    dedup_experimental(&mut paths, 150);
    assert_eq!(paths.len(), 1);
    let result = aligner.arbitrate(&query, paths, 30, true).unwrap();
    assert_eq!(crate::mapq::mapq_of(&result), 60);
}

#[test]
fn exact_paths_disabled_restores_original_annotation_preference() {
    let (index, query, library, pos, cigar) = terminal_fixture(2, false, false);
    let mut config = AlignConfig::rna_default();
    config.exact_path_competition = false;
    let mut aligner = Aligner::new(&index, config).with_lib(library);
    let mut continuous = placement(200, vec![CigarOp::Match(150)], Strand::Plus);
    let mut spliced = placement(pos, cigar.clone(), Strand::Plus);
    aligner.rescore_strict(&query, &mut continuous);
    aligner.rescore_strict(&query, &mut spliced);
    let mut paths = vec![continuous, spliced];
    dedup_by_footprint_with_options(&mut paths, 150, false);
    assert_eq!(paths.len(), 1);
    let result = aligner.arbitrate(&query, paths, 30, true).unwrap();
    assert_eq!(result.cigar, cigar);
    assert_eq!(result.score, 306);
    assert_eq!(crate::mapq::mapq_of(&result), 60);
}

pub(crate) fn native_l1_fixture() -> (
    Index<'static>,
    Vec<u8>,
    Arc<JunctionLib>,
    Arc<esperanto_engine::L1Index>,
) {
    let (index, query, library, _, _) = terminal_fixture(2, false, false);
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let serial = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let directory =
        std::env::temp_dir().join(format!("esp-exact-path-l1-{}-{serial}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut genome = Vec::new();
    index.reference.contigs[0].decode_append(0, 2000, &mut genome);
    let fasta = directory.join("reference.fa");
    std::fs::write(
        &fasta,
        [b">chrT\n".to_vec(), genome, b"\n".to_vec()].concat(),
    )
    .unwrap();
    std::fs::write(
        directory.join("reference.fa.fai"),
        b"chrT\t2000\t6\t2000\t2001\n",
    )
    .unwrap();
    let gtf = directory.join("annotation.gtf");
    std::fs::write(
        &gtf,
        concat!(
            "chrT\tx\ttranscript\t151\t450\t.\t+\t.\ttranscript_id \"a\";\n",
            "chrT\tx\texon\t151\t450\t.\t+\t.\ttranscript_id \"a\";\n",
            "chrT\tx\ttranscript\t151\t1200\t.\t+\t.\ttranscript_id \"b\";\n",
            "chrT\tx\texon\t151\t348\t.\t+\t.\ttranscript_id \"b\";\n",
            "chrT\tx\texon\t1101\t1200\t.\t+\t.\ttranscript_id \"b\";\n",
        ),
    )
    .unwrap();
    let tidx = directory.join("reference.tidx");
    esperanto_tidx::build(
        &gtf,
        &fasta,
        &tidx,
        &esperanto_tidx::BuildOptions::default(),
    )
    .unwrap();
    let l1 = Arc::new(esperanto_engine::L1Index::build(&tidx, &gtf, &fasta).unwrap());
    (index, query, library, l1)
}

#[test]
fn exact_paths_native_l1_cache_and_mate_preserve_alternatives() {
    let (index, query, library, l1) = native_l1_fixture();
    for strand in [Strand::Plus, Strand::Minus] {
        let raw = if strand == Strand::Plus {
            query.clone()
        } else {
            revcomp(&query)
        };
        let mut aligner = Aligner::new(&index, experimental_config())
            .with_lib(library.clone())
            .with_l1(l1.clone());
        let (result, cache) = aligner.align_read_full(&raw);
        let result = result.unwrap();
        assert_eq!(result.cigar, vec![CigarOp::Match(150)]);
        assert_eq!(crate::mapq::mapq_of(&result), 0);
        assert!(result.exact_path_ambiguous);
        let cache = cache.unwrap();
        let pool = cache.diagnostic_pool(&index, 150, "test");
        assert!(pool.candidates.iter().any(|p| p.cigar == "148M752N2M"));
        assert!(pool.candidates.iter().any(|p| p.cigar == "150M"));
        let mut rescued = aligner
            .align_read_with_mate(&raw, 0, 550, Some(&cache))
            .unwrap();
        assert!(rescued.exact_path_ambiguous);
        rescued.mapq_override = Some(60);
        assert_eq!(crate::mapq::mapq_of(&rescued), 0);
    }
}

#[test]
fn exact_paths_unknown_bases_and_partial_paths_do_not_claim_full_exactness() {
    let (index, mut query, library, pos, cigar) = terminal_fixture(2, false, false);
    query[50] = b'N';
    let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library);
    let mut a = placement(200, vec![CigarOp::Match(150)], Strand::Plus);
    let mut b = placement(pos, cigar, Strand::Plus);
    aligner.rescore_strict(&query, &mut a);
    aligner.rescore_strict(&query, &mut b);
    let mut paths = vec![a, b];
    dedup_experimental(&mut paths, 150);
    assert_eq!(paths.len(), 1);
    let result = aligner.arbitrate(&query, paths, 30, true).unwrap();
    assert!(!result.exact_path_ambiguous);

    query.clear();
    index.reference.contigs[0].decode_append(200, 350, &mut query);
    let mut a = placement(200, vec![CigarOp::Match(150)], Strand::Plus);
    let mut b = placement(
        200,
        vec![CigarOp::Match(149), CigarOp::SoftClip(1)],
        Strand::Plus,
    );
    aligner.rescore_strict(&query, &mut a);
    aligner.rescore_strict(&query, &mut b);
    let mut paths = vec![a, b];
    dedup_experimental(&mut paths, 150);
    assert_eq!(paths.len(), 1);
    let result = aligner.arbitrate(&query, paths, 30, true).unwrap();
    assert!(!result.exact_path_ambiguous);
}

#[test]
fn exact_paths_mate_coordinate_window_cannot_erase_overlapping_splice_alternative() {
    let mut genome = random(5000);
    let query = genome[3000..3150].to_vec();
    genome[100..102].copy_from_slice(&query[..2]);
    let fasta = [b">chrT\n".to_vec(), genome, b"\n".to_vec()].concat();
    let index = Index::build(
        parse_fasta_bytes(&fasta).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    let library = Arc::new(JunctionLib::build([Junction {
        contig: 0,
        start: 102,
        end: 3002,
        minus_strand: false,
    }]));
    let mut aligner = Aligner::new(&index, experimental_config()).with_lib(library);
    let mut a = placement(3000, vec![CigarOp::Match(150)], Strand::Plus);
    let mut b = placement(
        100,
        vec![
            CigarOp::Match(2),
            CigarOp::RefSkip(2900),
            CigarOp::Match(148),
        ],
        Strand::Plus,
    );
    aligner.rescore_strict(&query, &mut a);
    aligner.rescore_strict(&query, &mut b);
    let mut paths = vec![a, b];
    dedup_experimental(&mut paths, 150);
    let cache = CachedPlacements {
        placements: paths,
        n_seeds: 30,
        genomic_any: true,
        seed_budget: None,
    };
    let result = aligner
        .align_read_with_mate(&query, 0, 3200, Some(&cache))
        .unwrap();
    assert!(
        result.exact_path_ambiguous,
        "overlapping exact path is not an off-locus competitor"
    );
    assert_eq!(crate::mapq::mapq_of(&result), 0);
}

#[test]
fn default_policy_preserves_annotation_arbitration_in_native_l1_and_cache() {
    let (index, query, library, l1) = native_l1_fixture();
    for strand in [Strand::Plus, Strand::Minus] {
        let raw = if strand == Strand::Plus {
            query.clone()
        } else {
            revcomp(&query)
        };
        let mut aligner = Aligner::new(&index, AlignConfig::rna_default())
            .with_lib(library.clone())
            .with_l1(l1.clone());
        let (result, cache) = aligner.align_read_full(&raw);
        let result = result.unwrap();
        assert_eq!(
            result.cigar,
            vec![
                CigarOp::Match(148),
                CigarOp::RefSkip(752),
                CigarOp::Match(2)
            ],
            "default must restore pre-experiment annotation arbitration"
        );
        assert!(!result.exact_path_ambiguous);
        assert_eq!(crate::mapq::mapq_of(&result), 60);
        let result = aligner
            .align_read_with_mate(&raw, 0, 550, cache.as_ref())
            .unwrap();
        assert!(!result.exact_path_ambiguous);
        assert_eq!(
            result.cigar,
            vec![
                CigarOp::Match(148),
                CigarOp::RefSkip(752),
                CigarOp::Match(2)
            ]
        );
    }
}
