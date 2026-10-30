//! Real paired pipeline controls for endpoint/splice and two-pass interactions.
use super::*;

#[test]
fn pair_substitution_native_sparse_terminal_edit_competes_with_clipping() {
    let mut genome = crate::fragment::tests::body(8000);
    let mut body = genome[2000..2131].to_vec();
    // Both sequenced mates have Phred>=20 at the three internal positions
    // under this helper's deliberately nonuniform quality pattern.
    for position in [1, 21, 61, 62] {
        genome[2000 + position] = b'T';
        body[position] = b'C';
    }
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let pair = crate::fragment::tests::pair(&body, [150, 150], 0);
    let result = run_with_paths(
        "pair-substitution-red",
        &index,
        &[pair],
        true,
        None,
        false,
        true,
        true,
        false,
        None,
    );
    for r in &result.records {
        let aligned: usize = r
            .cigar()
            .iter()
            .map(|op| {
                let op = op.unwrap();
                if op.kind().consumes_read() && op.kind().consumes_reference() {
                    op.len()
                } else {
                    0
                }
            })
            .sum();
        assert_eq!(
            aligned, 131,
            "corroborated sparse terminal substitution must compete in the pair-supported domain"
        );
    }
}

#[test]
fn pair_substitution_native_directions_orders_passes_and_ablation() {
    for ag in [false, true] {
        for terminal in [1, 129] {
            let mut genome = crate::fragment::tests::body(8000);
            let mut body = genome[2000..2131].to_vec();
            for q in [terminal, 21, 61, 62] {
                genome[2000 + q] = if ag { b'A' } else { b'T' };
                body[q] = if ag { b'G' } else { b'C' };
            }
            let reference =
                crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat())
                    .unwrap();
            let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
            for swapped in [false, true] {
                // Reverse the molecule while retaining each sequencing mate's
                // own adapter family; swapping entire FASTQs breaks that evidence.
                let fragment = if swapped {
                    crate::align::revcomp(&body)
                } else {
                    body.clone()
                };
                let pair = crate::fragment::tests::pair(&fragment, [150, 150], 0);
                for two_pass in [false, true] {
                    for enabled in [false, true] {
                        let result = run_with_pair_substitution(
                            &format!("pair-subst-{ag}-{terminal}-{swapped}-{two_pass}-{enabled}"),
                            &index,
                            &[pair.clone()],
                            true,
                            None,
                            two_pass,
                            true,
                            true,
                            false,
                            None,
                            true,
                            enabled,
                            Some(41),
                        );
                        for (mate, r) in result.records.iter().enumerate() {
                            let aligned: usize = r
                                .cigar()
                                .iter()
                                .map(|op| {
                                    let op = op.unwrap();
                                    if op.kind().consumes_read() && op.kind().consumes_reference() {
                                        op.len()
                                    } else {
                                        0
                                    }
                                })
                                .sum();
                            assert_eq!(aligned, if enabled { 131 } else { 129 },
                                "ag={ag} terminal={terminal} swapped={swapped} pass2={two_pass} enabled={enabled}");
                            let reverse = r.flags().is_reverse_complemented();
                            let expected = if reverse {
                                crate::align::revcomp(&pair[mate])
                            } else {
                                pair[mate].clone()
                            };
                            assert_eq!(r.sequence().iter().collect::<Vec<_>>(), expected);
                            assert_eq!(
                                r.quality_scores().iter().collect::<Vec<_>>(),
                                vec![41; 150]
                            );
                            if enabled {
                                assert_eq!(r.mapping_quality().unwrap().get(), 0);
                                assert_eq!(r.alignment_start().unwrap().unwrap().get(), 2001);
                                assert!(format!("{:?}", r.data().get(b"RE"))
                                    .contains("pair_substitution"));
                                assert_eq!(
                                    result.rows[mate]["mate_attempt"]["pair_substitution"]
                                        ["decision"],
                                    "extended_mapq0"
                                );
                            }
                        }
                        assert!(result.discoveries.is_empty());
                        assert!(result.library.is_empty());
                    }
                }
            }
        }
    }
}

#[test]
fn pair_substitution_native_quality_endpoint_and_cache_controls() {
    let mut genome = crate::fragment::tests::body(8000);
    let mut body = genome[2000..2131].to_vec();
    for q in [1, 21, 61, 62] {
        genome[2000 + q] = b'T';
        body[q] = b'C';
    }
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let pair = crate::fragment::tests::pair(&body, [150, 150], 0);
    let pairs = [
        pair.clone(),
        [pair[0].clone(), vec![b'N'; 150]],
        pair.clone(),
    ];
    for two_pass in [false, true] {
        let result = run_with_pair_substitution(
            &format!("pair-subst-cache-{two_pass}"),
            &index,
            &pairs,
            true,
            None,
            two_pass,
            true,
            true,
            false,
            None,
            true,
            true,
            Some(41),
        );
        for i in [0, 2] {
            assert!(format!("{:?}", result.records[2 * i].data().get(b"RE"))
                .contains("pair_substitution"));
        }
        assert!(!format!("{:?}", result.records[2].data().get(b"RE")).contains("pair_substitution"));
        for (endpoints, quality) in [(true, 19), (false, 41)] {
            let result = run_with_pair_substitution(
                &format!("pair-subst-no-witness-{two_pass}-{endpoints}-{quality}"),
                &index,
                &[pair.clone()],
                endpoints,
                None,
                two_pass,
                true,
                true,
                false,
                None,
                true,
                true,
                Some(quality),
            );
            for r in &result.records {
                assert!(!format!("{:?}", r.data().get(b"RE")).contains("pair_substitution"));
            }
        }
    }
}

#[test]
fn boundary_core_native_frames_orders_threads_and_ablation() {
    let mut genome = crate::fragment::tests::body(8000);
    genome[2076..2080].fill(b'A');
    genome[2075] = b'C';
    // A concrete negative-score flank prevents a different mismatching DP
    // extension from replacing the exact source path this test exercises.
    genome[2080..2230].fill(b'C');
    let query = [&genome[2000..2080], &vec![b'A'; 70]].concat();
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap(),
        crate::seed::SeedParams { k: 15, w: 5 },
    );
    for minus in [false, true] {
        let raw = if minus {
            crate::align::revcomp(&query)
        } else {
            query.clone()
        };
        let mate = if minus {
            genome[2300..2450].to_vec()
        } else {
            crate::align::revcomp(&genome[2300..2450])
        };
        for swapped in [false, true] {
            let pair = if swapped {
                [mate.clone(), raw.clone()]
            } else {
                [raw.clone(), mate.clone()]
            };
            let target = usize::from(swapped);
            for two_pass in [false, true] {
                for threads in [1, 3] {
                    for enabled in [false, true] {
                        let out = run_with_boundary_core(
                            &format!(
                                "core-native-{minus}-{swapped}-{two_pass}-{threads}-{enabled}"
                            ),
                            &index,
                            &[pair.clone()],
                            false,
                            None,
                            two_pass,
                            true,
                            true,
                            false,
                            None,
                            true,
                            true,
                            None,
                            enabled,
                            threads,
                        );
                        let r = &out.records[target];
                        assert_eq!(!r.flags().is_unmapped(), enabled,
                            "minus={minus} swapped={swapped} two_pass={two_pass} threads={threads} enabled={enabled} rows={:?}",
                            out.rows);
                        if enabled {
                            assert_eq!(r.alignment_start().unwrap().unwrap().get(), 2001);
                            assert_eq!(r.mapping_quality().unwrap().get(), 0);
                            assert_eq!(r.flags().is_reverse_complemented(), minus);
                            assert_eq!(r.sequence().iter().collect::<Vec<_>>(), query);
                            assert!(format!("{:?}", r.data().get(b"XB").unwrap().unwrap())
                                .contains("query_start=76;length=4;assignment=unknown"));
                            let mut qualities: Vec<_> =
                                (0..150).map(|q| ((q + target) % 41) as u8).collect();
                            if minus {
                                qualities.reverse();
                            }
                            assert_eq!(r.quality_scores().iter().collect::<Vec<_>>(), qualities);
                            assert!(r.flags().is_properly_segmented());
                            assert!(!out.rows[target]["mate_attempt"]["boundary_core_proposals"]
                                .as_array()
                                .unwrap()
                                .is_empty());
                        } else {
                            assert!(r.data().get(b"XB").is_none());
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn boundary_core_fully_biological_mapped_run_is_preserved() {
    let mut genome = crate::fragment::tests::body(8000);
    genome[2080..2150].fill(b'A');
    genome[2079] = b'C';
    genome[2150] = b'C';
    let query = genome[2000..2150].to_vec();
    let mate = crate::align::revcomp(&genome[2300..2450]);
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap(),
        crate::seed::SeedParams { k: 15, w: 5 },
    );
    for enabled in [false, true] {
        let out = run_with_boundary_core(
            &format!("core-biological-{enabled}"),
            &index,
            &[[query.clone(), mate.clone()]],
            false,
            None,
            true,
            true,
            true,
            false,
            None,
            true,
            true,
            None,
            enabled,
            2,
        );
        let r = &out.records[0];
        assert!(!r.flags().is_unmapped());
        assert_eq!(r.sequence().iter().collect::<Vec<_>>(), query);
        assert_eq!(r.alignment_start().unwrap().unwrap().get(), 2001);
        assert_eq!(
            r.cigar().iter().map(|op| op.unwrap().len()).sum::<usize>(),
            150
        );
        assert!(r
            .cigar()
            .iter()
            .all(|op| op.unwrap().kind().consumes_reference()));
        assert!(r.data().get(b"XB").is_none());
    }
}

#[test]
fn short_body_conditional_native_remote_repeat_mate_orders_and_passes() {
    let mut genome = crate::fragment::tests::body(8000);
    let body = genome[2000..2080].to_vec();
    genome[5000..5080].copy_from_slice(&body);
    let tail = *b"ACGT"
        .iter()
        .find(|&&b| b != body[79] && b != genome[2080])
        .unwrap();
    let query = [body.as_slice(), &vec![tail; 70]].concat();
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    for minus in [false, true] {
        let raw = if minus {
            crate::align::revcomp(&query)
        } else {
            query.clone()
        };
        let mate = if minus {
            genome[2300..2450].to_vec()
        } else {
            crate::align::revcomp(&genome[2300..2450])
        };
        for swapped in [false, true] {
            let pair = if swapped {
                [mate.clone(), raw.clone()]
            } else {
                [raw.clone(), mate.clone()]
            };
            let target = usize::from(swapped);
            for two_pass in [false, true] {
                let observed = run_with_short_body(
                    &format!("conditional-native-{minus}-{swapped}-{two_pass}"),
                    &index,
                    &[pair.clone()],
                    false,
                    None,
                    two_pass,
                    true,
                    true,
                    false,
                    None,
                    true,
                );
                let record = &observed.records[target];
                assert!(!record.flags().is_unmapped());
                assert_eq!(record.alignment_start().unwrap().unwrap().get(), 2001);
                assert_eq!(record.mapping_quality().unwrap().get(), 0);
                assert_eq!(record.flags().is_reverse_complemented(), minus);
                assert_eq!(record.sequence().iter().collect::<Vec<_>>(), query);
                assert_eq!(
                    observed.rows[target]["mate_attempt"]["short_body"]["supported_candidates"][0]
                        ["support_scope"],
                    "mate_conditioned_index",
                );
            }
        }
    }
}

#[test]
fn short_body_native_pair_recovers_visible_boundary_without_changing_raw_reads() {
    let genome = crate::fragment::tests::body(8000);
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">chrT\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let body = &genome[2000..2080];
    let tail = *b"ACGT"
        .iter()
        .find(|&&b| b != body[79] && b != genome[2080])
        .unwrap();
    let query = [body, &vec![tail; 70]].concat();
    for minus in [false, true] {
        let raw = if minus {
            crate::align::revcomp(&query)
        } else {
            query.clone()
        };
        let mate = if minus {
            genome[2300..2450].to_vec()
        } else {
            crate::align::revcomp(&genome[2300..2450])
        };
        for swapped in [false, true] {
            let pair = if swapped {
                [mate.clone(), raw.clone()]
            } else {
                [raw.clone(), mate.clone()]
            };
            let target = usize::from(swapped);
            for two_pass in [false, true] {
                for enabled in [false, true] {
                    let observed = run_with_short_body(
                        &format!("shortbody-{minus}-{swapped}-{two_pass}-{enabled}"),
                        &index,
                        &[pair.clone()],
                        false,
                        None,
                        two_pass,
                        true,
                        true,
                        false,
                        None,
                        enabled,
                    );
                    let record = &observed.records[target];
                    assert_eq!(
                        !record.flags().is_unmapped(), enabled,
                        "exact 80bp body with visible 70bp tail and unique proper mate must qualify"
                    );
                    let sequence = if enabled { query.clone() } else { raw.clone() };
                    assert_eq!(record.sequence().iter().collect::<Vec<_>>(), sequence);
                    let mut qualities: Vec<_> =
                        (0..150).map(|q| ((q + target) % 41) as u8).collect();
                    if enabled && minus {
                        qualities.reverse();
                    }
                    assert_eq!(
                        record.quality_scores().iter().collect::<Vec<_>>(),
                        qualities
                    );
                    if enabled {
                        assert_eq!(record.alignment_start().unwrap().unwrap().get(), 2001);
                        assert_eq!(record.mapping_quality().unwrap().get(), 0);
                        assert!(record.flags().is_properly_segmented());
                        assert_eq!(record.flags().is_reverse_complemented(), minus);
                        assert_eq!(
                            observed.rows[target]["mate_attempt"]["short_body"]["decision"],
                            "qualified_mapq0"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn short_body_native_cache_is_pair_specific_and_full_genomic_tail_stays_mapped() {
    let mut genome = crate::fragment::tests::body(8000);
    let body = genome[2000..2080].to_vec();
    let tail = *b"ACGT"
        .iter()
        .find(|&&b| b != body[79] && b != genome[2080])
        .unwrap();
    let query = [body.as_slice(), &vec![tail; 70]].concat();
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let good = crate::align::revcomp(&genome[2300..2450]);
    let pairs = [
        [query.clone(), good.clone()],
        [query.clone(), genome[2300..2450].to_vec()], // wrong same strand
        [query.clone(), vec![b'N'; 150]],             // no mate
        [query.clone(), crate::align::revcomp(&genome[6000..6150])], // excessive span
        [query.clone(), good],                        // same cache key, valid pair again
    ];
    for two_pass in [false, true] {
        let observed = run_with_short_body(
            &format!("shortbody-cache-{two_pass}"),
            &index,
            &pairs,
            false,
            None,
            two_pass,
            true,
            true,
            false,
            None,
            true,
        );
        for i in 0..pairs.len() {
            assert_eq!(
                !observed.records[2 * i].flags().is_unmapped(),
                i == 0 || i == 4
            );
        }
    }
    // A true genomic homopolymer remains a normal full alignment. The fallback
    // cannot replace an already mapped result with a truncated interpretation.
    genome[2080..2150].fill(tail);
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">c\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    for endpoints in [false, true] {
        let pair = [query.clone(), crate::align::revcomp(&genome[2300..2450])];
        let on = run_with_short_body(
            &format!("shortbody-full-on-{endpoints}"),
            &index,
            &[pair.clone()],
            endpoints,
            None,
            false,
            true,
            true,
            false,
            None,
            true,
        );
        let off = run_with_short_body(
            &format!("shortbody-full-off-{endpoints}"),
            &index,
            &[pair],
            endpoints,
            None,
            false,
            true,
            true,
            false,
            None,
            false,
        );
        assert!(!on.records[0].flags().is_unmapped());
        assert_eq!(
            on.records[0].cigar().as_ref(),
            off.records[0].cigar().as_ref()
        );
        assert_eq!(
            on.records[0].mapping_quality(),
            off.records[0].mapping_quality()
        );
        assert_eq!(on.records[0].sequence().iter().collect::<Vec<_>>(), query);
    }
}

#[test]
fn short_body_native_splice_final_pass_does_not_change_discovery_library() {
    let (index, body, junction) = splice_world(80);
    let tail = *b"ACGT".iter().find(|&&b| b != body[79]).unwrap();
    let query = [body.as_slice(), &vec![tail; 70]].concat();
    let mate = crate::align::revcomp(&index.reference.contigs[0].slice_ascii(9200, 9350));
    let lib = Arc::new(JunctionLib::build([junction]));
    let pair = [query.clone(), mate];
    for two_pass in [false, true] {
        let on = run_with_short_body(
            &format!("shortbody-splice-on-{two_pass}"),
            &index,
            &[pair.clone(), pair.clone()],
            false,
            Some(lib.clone()),
            two_pass,
            true,
            true,
            false,
            None,
            true,
        );
        let off = run_with_short_body(
            &format!("shortbody-splice-off-{two_pass}"),
            &index,
            &[pair.clone(), pair.clone()],
            false,
            Some(lib.clone()),
            two_pass,
            true,
            true,
            false,
            None,
            false,
        );
        if two_pass {
            assert_eq!(on.library, off.library);
            assert_eq!(
                format!("{:?}", on.discoveries),
                format!("{:?}", off.discoveries)
            );
        }
        assert!(!on.records[0].flags().is_unmapped());
        assert!(off.records[0].flags().is_unmapped());
        assert_eq!(on.records[0].mapping_quality().unwrap().get(), 0);
        assert!(on.rows[0]["mate_attempt"]["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .any(|p| p["cigar"] == "40M1000N40M70S"));
    }
}

#[derive(Clone)]
struct Sink(Arc<std::sync::Mutex<Vec<u8>>>);
impl Write for Sink {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

struct Observed {
    records: Vec<noodles::bam::Record>,
    rows: Vec<serde_json::Value>,
    discoveries: Vec<FilteredDiscovery>,
    library: Vec<(Junction, u32)>,
}

fn run(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
) -> Observed {
    run_with_context(name, index, pairs, enabled, lib, two_pass, true)
}

#[allow(clippy::too_many_arguments)]
fn run_with_context(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
) -> Observed {
    run_with_background(
        name,
        index,
        pairs,
        enabled,
        lib,
        two_pass,
        repeat_context,
        true,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_background(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
    context_background: bool,
) -> Observed {
    run_with_paths(
        name,
        index,
        pairs,
        enabled,
        lib,
        two_pass,
        repeat_context,
        context_background,
        true,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_paths(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
    context_background: bool,
    exact_path_competition: bool,
    l1: Option<Arc<esperanto_engine::L1Index>>,
) -> Observed {
    run_with_short_body(
        name,
        index,
        pairs,
        enabled,
        lib,
        two_pass,
        repeat_context,
        context_background,
        exact_path_competition,
        l1,
        true,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_short_body(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
    context_background: bool,
    exact_path_competition: bool,
    l1: Option<Arc<esperanto_engine::L1Index>>,
    mate_short_body: bool,
) -> Observed {
    run_with_pair_substitution(
        name,
        index,
        pairs,
        enabled,
        lib,
        two_pass,
        repeat_context,
        context_background,
        exact_path_competition,
        l1,
        mate_short_body,
        true,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_pair_substitution(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
    context_background: bool,
    exact_path_competition: bool,
    l1: Option<Arc<esperanto_engine::L1Index>>,
    mate_short_body: bool,
    pair_substitution: bool,
    phred: Option<u8>,
) -> Observed {
    run_with_boundary_core(
        name,
        index,
        pairs,
        enabled,
        lib,
        two_pass,
        repeat_context,
        context_background,
        exact_path_competition,
        l1,
        mate_short_body,
        pair_substitution,
        phred,
        true,
        2,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_boundary_core(
    name: &str,
    index: &Index,
    pairs: &[[Vec<u8>; 2]],
    enabled: bool,
    lib: Option<Arc<JunctionLib>>,
    two_pass: bool,
    repeat_context: bool,
    context_background: bool,
    exact_path_competition: bool,
    l1: Option<Arc<esperanto_engine::L1Index>>,
    mate_short_body: bool,
    pair_substitution: bool,
    phred: Option<u8>,
    boundary_core: bool,
    threads: usize,
) -> Observed {
    let dir = std::env::temp_dir().join(format!("esp-fragment-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let paths = [dir.join("r1.fq"), dir.join("r2.fq")];
    for mate in 0..2 {
        let mut text = Vec::new();
        for (i, pair) in pairs.iter().enumerate() {
            let seq = &pair[mate];
            text.extend_from_slice(format!("@p{i}/{}\n", mate + 1).as_bytes());
            text.extend_from_slice(seq);
            text.extend_from_slice(b"\n+\n");
            text.extend((0..seq.len()).map(|q| b'!' + phred.unwrap_or(((q + mate) % 41) as u8)));
            text.push(b'\n');
        }
        std::fs::write(&paths[mate], text).unwrap();
    }
    let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
    let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut config = AlignConfig::rna_default();
    config.fragment_endpoints = enabled;
    config.repeat_context = repeat_context;
    config.context_background = context_background;
    config.exact_path_competition = exact_path_competition;
    config.mate_short_body = mate_short_body;
    config.pair_substitution = pair_substitution;
    config.boundary_core = boundary_core;
    let mut out = PipelineOut {
        bam: Some(Box::new(Sink(raw.clone()))),
        unmapped_fq: Box::new(std::io::sink()),
        diagnostics: Some(Box::new(Sink(diag.clone()))),
        index,
        config,
        jlib: lib,
        l1,
    };
    let (discoveries, library) = if two_pass {
        let (_, kept, library, _) = run_pe_2pass(&mut out, &paths[0], &paths[1], threads).unwrap();
        (kept, library)
    } else {
        run_pe(&mut out, &paths[0], &paths[1], threads).unwrap();
        (Vec::new(), Vec::new())
    };
    drop(out);
    let bytes = raw.lock().unwrap().clone();
    let mut reader = noodles::bam::io::Reader::new(std::io::Cursor::new(bytes));
    reader.read_header().unwrap();
    let records = reader.records().map(Result::unwrap).collect();
    let rows = std::str::from_utf8(&diag.lock().unwrap())
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    Observed {
        records,
        rows,
        discoveries,
        library,
    }
}

#[test]
fn exact_paths_native_paired_bam_keeps_ambiguity_and_raw_frame() {
    let (index, query, lib, l1) = crate::align::exact_path_tests::native_l1_fixture();
    let mut mate = Vec::new();
    index.reference.contigs[0].decode_append(550, 700, &mut mate);
    for minus in [false, true] {
        let pair = if minus {
            [crate::align::revcomp(&query), mate.clone()]
        } else {
            [query.clone(), crate::align::revcomp(&mate)]
        };
        for two_pass in [false, true] {
            for enabled in [false, true] {
                let result = run_with_paths(
                    &format!("exact-path-{minus}-{two_pass}-{enabled}"),
                    &index,
                    &[pair.clone()],
                    true,
                    Some(lib.clone()),
                    two_pass,
                    true,
                    true,
                    enabled,
                    Some(l1.clone()),
                );
                let record = &result.records[0];
                assert_eq!(
                    record.mapping_quality().unwrap().get(),
                    if enabled { 0 } else { 60 }
                );
                let evidence = format!("{:?}", record.data().get(b"EK").unwrap().unwrap());
                assert_eq!(evidence.contains("exact_path_tie"), enabled);
                let sequence = if minus {
                    crate::align::revcomp(&pair[0])
                } else {
                    pair[0].clone()
                };
                assert_eq!(record.sequence().iter().collect::<Vec<_>>(), sequence);
                let mut qualities: Vec<_> = (0..150).map(|q| (q % 41) as u8).collect();
                if minus {
                    qualities.reverse();
                }
                assert_eq!(
                    record.quality_scores().iter().collect::<Vec<_>>(),
                    qualities
                );
                assert_eq!(
                    result.rows[0]["independent"]["candidates"]
                        .as_array()
                        .unwrap()
                        .len(),
                    if enabled { 2 } else { 1 }
                );
            }
        }
    }
}

#[test]
fn intron_union_native_paired_splice_has_nonzero_signed_template_span() {
    let genome = crate::fragment::tests::body(4000);
    let fasta = [b">chrUnion\n".to_vec(), genome.clone(), b"\n".to_vec()].concat();
    let index = Index::build(
        crate::fasta::parse_fasta_bytes(&fasta).unwrap(),
        crate::seed::SeedParams { k: 15, w: 5 },
    );
    // Transcript[200,300)+[2300,2400), two overlapping150bp mates.
    let r1 = [genome[200..300].to_vec(), genome[2300..2350].to_vec()].concat();
    let r2 = [genome[250..300].to_vec(), genome[2300..2400].to_vec()].concat();
    let lib = Arc::new(JunctionLib::build([Junction {
        contig: 0,
        start: 300,
        end: 2300,
        minus_strand: false,
    }]));
    for reverse_order in [false, true] {
        let pair = if reverse_order {
            [crate::align::revcomp(&r2), r1.clone()]
        } else {
            [r1.clone(), crate::align::revcomp(&r2)]
        };
        for two_pass in [false, true] {
            let result = run_with_paths(
                &format!("intron-union-{reverse_order}-{two_pass}"),
                &index,
                &[pair.clone()],
                false,
                Some(lib.clone()),
                two_pass,
                true,
                true,
                false,
                None,
            );
            assert_eq!(result.records.len(), 2);
            for (mate, record) in result.records.iter().enumerate() {
                assert!(!record.flags().is_unmapped());
                assert!(record.flags().is_properly_segmented());
                assert_eq!(
                    record.template_length(),
                    if (mate == 0) == reverse_order {
                        -200
                    } else {
                        200
                    },
                    "shared intron must not produce zero template length"
                );
                let is_left = (mate == 0) != reverse_order;
                assert_eq!(
                    record.alignment_start().unwrap().unwrap().get(),
                    if is_left { 201 } else { 251 }
                );
                assert_eq!(
                    record.mate_alignment_start().unwrap().unwrap().get(),
                    if is_left { 251 } else { 201 }
                );
                let expected = if is_left { r1.clone() } else { r2.clone() };
                assert_eq!(record.sequence().iter().collect::<Vec<_>>(), expected);
                let mut qualities: Vec<_> = (0..150).map(|q| ((q + mate) % 41) as u8).collect();
                if !is_left {
                    qualities.reverse();
                }
                assert_eq!(
                    record.quality_scores().iter().collect::<Vec<_>>(),
                    qualities
                );
            }
        }
    }
}

#[test]
fn context_background_paired_refusal_preserves_unmapped_raw_evidence() {
    let (index, fragment) = crate::align::background_tests::fixture(false);
    let pair = crate::fragment::tests::pair(&fragment, [150, 150], 0);
    let off = run_with_background(
        "background-off",
        &index,
        &[pair.clone()],
        true,
        None,
        false,
        true,
        false,
    );
    assert!(off.records.iter().all(|rec| !rec.flags().is_unmapped()));
    assert!(off
        .rows
        .iter()
        .all(|row| row["independent"]["seed_budget"]["repeat_context"]
            .get("background")
            .is_none()));
    assert_eq!(
        crate::fragment::infer(&pair[0], &pair[1])
            .unwrap()
            .fragment_length,
        118
    );
    for two_pass in [false, true] {
        let name = if two_pass {
            "background-null-two-pass"
        } else {
            "background-null-one-pass"
        };
        let on = run_with_context(name, &index, &[pair.clone()], true, None, two_pass, true);
        assert_eq!((on.records.len(), on.rows.len()), (2, 2));
        for (mate, (rec, row)) in on.records.iter().zip(&on.rows).enumerate() {
            assert!(
                rec.flags().is_unmapped(),
                "background-dominated pair must stay unresolved"
            );
            assert_eq!(row["pass2_primary_placed"], false);
            assert_eq!(row["effective_read_length"], 118);
            let bg = &row["independent"]["seed_budget"]["repeat_context"]["background"];
            assert!(bg["score"].as_i64().unwrap() >= 152);
            assert!(bg["rejected_candidates"].as_u64().unwrap() > 0);
            assert_eq!(rec.sequence().iter().collect::<Vec<_>>(), pair[mate]);
            assert_eq!(
                rec.quality_scores().iter().collect::<Vec<_>>(),
                (0..150)
                    .map(|q| ((q + mate) % 41) as u8)
                    .collect::<Vec<_>>()
            );
            assert_eq!(rec.data().get(b"XF").unwrap().unwrap().as_int(), Some(118));
        }
    }
}

#[test]
fn ordinary_repeat_context_paired_bam_and_diagnostics_preserve_raw_frame() {
    let (index, fragment) = crate::align::repeat_context_tests::fixture(false, false, 6000);
    let pair = crate::fragment::tests::pair(&fragment, [150, 150], 0);
    assert_eq!(
        crate::fragment::infer(&pair[0], &pair[1])
            .unwrap()
            .fragment_length,
        112
    );
    let off = run_with_context(
        "repeat-context-off",
        &index,
        &[pair.clone()],
        true,
        None,
        false,
        false,
    );
    assert_eq!(off.records.len(), 2);
    assert_eq!(off.rows.len(), 2);
    assert!(off.records.iter().all(|r| r.flags().is_unmapped()));
    for two_pass in [false, true] {
        let name = if two_pass {
            "repeat-context-two-pass"
        } else {
            "repeat-context-one-pass"
        };
        let on = run_with_context(name, &index, &[pair.clone()], true, None, two_pass, true);
        assert_eq!(on.records.len(), 2);
        assert_eq!(on.rows.len(), 2);
        for (mate, (rec, row)) in on.records.iter().zip(&on.rows).enumerate() {
            assert!(!rec.flags().is_unmapped());
            assert_eq!(rec.reference_sequence_id().unwrap().unwrap(), 2);
            assert_eq!(rec.alignment_start().unwrap().unwrap().get(), 121);
            assert_eq!(rec.mapping_quality().unwrap().get(), 0);
            assert_eq!(row["mate"], mate + 1);
            assert_eq!(row["effective_read_length"], 112);
            assert_eq!(row["pass2_primary_placed"], true);
            assert_eq!(
                row["independent"]["seed_budget"]["repeat_context"]["query_length"],
                112
            );
            let ops: Vec<_> = rec
                .cigar()
                .iter()
                .map(|p| {
                    let p = p.unwrap();
                    (format!("{:?}", p.kind()), p.len())
                })
                .collect();
            assert_eq!(
                ops,
                if mate == 0 {
                    vec![("Match".to_owned(), 112), ("SoftClip".to_owned(), 38)]
                } else {
                    vec![("SoftClip".to_owned(), 38), ("Match".to_owned(), 112)]
                }
            );
            let mut seq: Vec<_> = rec.sequence().iter().collect();
            let mut qual: Vec<_> = rec.quality_scores().iter().collect();
            if rec.flags().is_reverse_complemented() {
                seq = crate::align::revcomp(&seq);
                qual.reverse();
            }
            assert_eq!(seq, pair[mate]);
            assert_eq!(
                qual,
                (0..150)
                    .map(|q| ((q + mate) % 41) as u8)
                    .collect::<Vec<_>>()
            );
            assert_eq!(rec.data().get(b"XF").unwrap().unwrap().as_int(), Some(112));
            assert_eq!(
                rec.data()
                    .iter()
                    .filter(|p| p.as_ref().unwrap().0 == *b"XF")
                    .count(),
                1
            );
        }
    }
}

#[test]
fn query_edge_completion_paired_pipeline_preserves_raw_frame() {
    let (index, fragment) = crate::align::prefix_tests::fixture();
    let pair = crate::fragment::tests::pair(&fragment, [150, 150], 0);
    assert_eq!(
        crate::fragment::infer(&pair[0], &pair[1])
            .unwrap()
            .fragment_length,
        142
    );
    let off = run(
        "edge-prefix-off",
        &index,
        &[pair.clone()],
        false,
        None,
        false,
    );
    assert_eq!(off.records.len(), 2);
    assert!(off.records.iter().all(|rec| !rec.flags().is_unmapped()));
    let on = run("edge-prefix-on", &index, &[pair.clone()], true, None, false);
    assert_eq!(on.records.len(), 2);
    assert_eq!(on.rows.len(), 2);
    for (mate, (rec, row)) in on.records.iter().zip(&on.rows).enumerate() {
        assert_eq!(row["mate"], mate + 1);
        assert!(
            !rec.flags().is_unmapped(),
            "effective prefix must retain its body anchor"
        );
        assert_eq!(rec.reference_sequence_id().unwrap().unwrap(), 2);
        assert_eq!(rec.alignment_start().unwrap().unwrap().get(), 121);
        assert_eq!(row["schema_version"], 2);
        assert_eq!(row["effective_read_length"], 142);
        assert_eq!(row["pass2_primary_placed"], true);
        let ops: Vec<_> = rec
            .cigar()
            .iter()
            .map(|op| {
                let op = op.unwrap();
                (format!("{:?}", op.kind()), op.len())
            })
            .collect();
        let expected = if mate == 0 {
            vec![("Match".to_owned(), 142), ("SoftClip".to_owned(), 8)]
        } else {
            vec![("SoftClip".to_owned(), 8), ("Match".to_owned(), 142)]
        };
        assert_eq!(ops, expected);
        let mut seq: Vec<_> = rec.sequence().iter().collect();
        let mut qual: Vec<_> = rec.quality_scores().iter().collect();
        if rec.flags().is_reverse_complemented() {
            seq = crate::align::revcomp(&seq);
            qual.reverse();
        }
        assert_eq!(seq, pair[mate]);
        assert_eq!(
            qual,
            (0..150)
                .map(|q| ((q + mate) % 41) as u8)
                .collect::<Vec<_>>()
        );
        assert_eq!(rec.data().get(b"XF").unwrap().unwrap().as_int(), Some(142));
        assert_eq!(
            rec.data()
                .iter()
                .filter(|field| field.as_ref().unwrap().0 == *b"XF")
                .count(),
            1
        );
    }
}

fn splice_world(n: usize) -> (Index, Vec<u8>, Junction) {
    // Independent random reference with a single canonical intron.
    let mut state = 0x72a984b1821u64;
    let mut genome: Vec<u8> = (0..20000)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect();
    genome[8000..8002].copy_from_slice(b"GT");
    genome[8998..9000].copy_from_slice(b"AG");
    let left = n / 2;
    let fragment = [
        genome[8000 - left..8000].to_vec(),
        genome[9000..9000 + n - left].to_vec(),
    ]
    .concat();
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">chrT\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let junction = Junction {
        contig: 0,
        start: 8000,
        end: 9000,
        minus_strand: false,
    };
    (index, fragment, junction)
}

#[test]
fn paired_fragment_splice_pipeline_preserves_both_mates_raw_frame() {
    let (index, fragment, junction) = splice_world(110);
    let pair = crate::fragment::tests::pair(&fragment, [150, 170], 0);
    let lib = Some(Arc::new(JunctionLib::build([junction])));
    let on = run(
        "splice-on",
        &index,
        &[pair.clone()],
        true,
        lib.clone(),
        false,
    );
    let off = run("splice-off", &index, &[pair.clone()], false, lib, false);
    assert_eq!(on.records.len(), 2);
    for (mate, (rec, row)) in on.records.iter().zip(&on.rows).enumerate() {
        assert!(!rec.flags().is_unmapped());
        assert_eq!(rec.alignment_start().unwrap().unwrap().get(), 7946);
        assert_eq!(row["schema_version"], 2);
        let expected = if mate == 0 {
            "55M1000N55M40S"
        } else {
            "60S55M1000N55M"
        };
        assert!(
            row["independent"]["candidates"]
                .as_array()
                .unwrap()
                .iter()
                .any(|c| c["cigar"] == expected && c["ordinary_mismatches"] == 0),
            "{row}"
        );
        let mut seq: Vec<u8> = rec.sequence().iter().collect();
        let mut qual: Vec<u8> = rec.quality_scores().iter().collect();
        if rec.flags().is_reverse_complemented() {
            seq = crate::align::revcomp(&seq);
            qual.reverse();
        }
        assert_eq!(seq, pair[mate]);
        assert_eq!(
            qual,
            (0..pair[mate].len())
                .map(|q| ((q + mate) % 41) as u8)
                .collect::<Vec<_>>()
        );
        assert!(format!("{:?}", rec.data()).contains("XF"));
    }
    assert!(off.rows.iter().all(|row| row["schema_version"] == 1));
}

#[test]
fn paired_fragment_two_pass_discovers_junction_without_tagging_unsupported_recipient() {
    let (index, fragment, junction) = splice_world(110);
    let pair = crate::fragment::tests::pair(&fragment, [200, 200], 0);
    // Raw110/200 coverage would fail. Two copies generate four confident
    // endpoint-supported observations; a separate untrimmed pair stays schema1.
    let contig = &index.reference.contigs[0];
    let recipient = [
        contig.slice_ascii(7898, 8000),
        contig.slice_ascii(9000, 9008),
    ]
    .concat();
    let recipient = [
        recipient,
        crate::align::revcomp(&contig.slice_ascii(9500, 9610)),
    ];
    let pairs = [pair.clone(), pair, recipient];
    let on = run("discovery-on", &index, &pairs, true, None, true);
    let off = run("discovery-off", &index, &pairs, false, None, true);
    assert!(
        on.discoveries
            .iter()
            .any(|(j, count, _)| *j == junction && *count >= 4),
        "{:?}",
        on.discoveries
    );
    assert!(on
        .library
        .iter()
        .any(|(j, count)| *j == junction && *count >= 4));
    let on_support = on.library.iter().find(|(j, _)| *j == junction).unwrap().1;
    let off_support = off
        .library
        .iter()
        .find(|(j, _)| *j == junction)
        .map_or(0, |x| x.1);
    assert!(on_support > off_support);
    assert_eq!(
        on.rows.len(),
        6,
        "pass1 must not emit duplicate diagnostics"
    );
    for (i, row) in on.rows.iter().enumerate() {
        assert_eq!(row["schema_version"], if i < 4 { 2 } else { 1 });
        assert!(row["pass2_primary_placed"].as_bool().unwrap());
    }
    let has_full_splice = |row: &serde_json::Value| {
        row["independent"]["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| {
                c["informative_query_bases"] == 110
                    && c["cigar"].as_str().unwrap().contains("1000N")
            })
    };
    assert!(
        has_full_splice(&on.rows[4]),
        "known junction must recover the unseeded short exon; on={} off={}",
        on.rows[4],
        off.rows[4]
    );
    let full_splice_score = |row: &serde_json::Value| {
        row["independent"]["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|c| {
                c["informative_query_bases"] == 110
                    && c["cigar"].as_str().unwrap().contains("1000N")
            })
            .map(|c| c["common_score"].as_i64().unwrap())
            .max()
    };
    assert_ne!(
        full_splice_score(&on.rows[4]),
        full_splice_score(&off.rows[4]),
        "new discovery must affect the unsupported recipient's candidate evidence"
    );
    // Explicit confidence gate: repeated MAPQ0 supporters cannot establish it.
    let (kept, _) = filter_discoveries(vec![(junction, 4, 0)], None);
    assert!(kept.is_empty());
}

#[test]
fn paired_fragment_long_dp_retains_terminal_edit_evidence_in_both_mates() {
    let mut reference_body = crate::fragment::tests::body(300);
    // The supported effective prefix exceeds256bp and exercises DP directly;
    // this control does not pretend to fix the separate short Myers tie.
    let positions: Vec<usize> = [0, 299].into_iter().chain((20..290).step_by(16)).collect();
    for &p in &positions {
        reference_body[p] = b'A';
    }
    let mut edited = reference_body.clone();
    for &p in &positions {
        edited[p] = b'G';
    }
    let genome = [vec![b'N'; 1000], reference_body, vec![b'N'; 1000]].concat();
    let reference =
        crate::fasta::parse_fasta_bytes(&[b">chrT\n".as_slice(), &genome, b"\n"].concat()).unwrap();
    let index = Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });
    let pair = crate::fragment::tests::pair(&edited, [500, 520], 0);
    let observed = run("terminal-dp", &index, &[pair.clone()], true, None, false);
    for (mate, row) in observed.rows.iter().enumerate() {
        let expected = if mate == 0 { "300M200S" } else { "220S300M" };
        assert!(row["pass2_primary_placed"].as_bool().unwrap(), "{row}");
        assert!(
            row["independent"]["candidates"]
                .as_array()
                .unwrap()
                .iter()
                .any(|c| c["cigar"] == expected
                    && c["edit_ag"] == positions.len()
                    && c["ordinary_mismatches"] == 0
                    && c["editing_score_mode"] == "ag"),
            "{row}"
        );
        assert!(!observed.records[mate].flags().is_unmapped());
    }
}
