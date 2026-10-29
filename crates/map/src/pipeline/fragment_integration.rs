//! Real paired pipeline controls for endpoint/splice and two-pass interactions.
use super::*;

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
            text.extend((0..seq.len()).map(|q| b'!' + ((q + mate) % 41) as u8));
            text.push(b'\n');
        }
        std::fs::write(&paths[mate], text).unwrap();
    }
    let raw = Arc::new(std::sync::Mutex::new(Vec::new()));
    let diag = Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut config = AlignConfig::rna_default();
    config.fragment_endpoints = enabled;
    let mut out = PipelineOut {
        bam: Some(Box::new(Sink(raw.clone()))),
        unmapped_fq: Box::new(std::io::sink()),
        diagnostics: Some(Box::new(Sink(diag.clone()))),
        index,
        config,
        jlib: lib,
        l1: None,
    };
    let (discoveries, library) = if two_pass {
        let (_, kept, library, _) = run_pe_2pass(&mut out, &paths[0], &paths[1], 2).unwrap();
        (kept, library)
    } else {
        run_pe(&mut out, &paths[0], &paths[1], 2).unwrap();
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
    assert!(off.records.iter().all(|rec| !rec.flags().is_unmapped()));
    let on = run("edge-prefix-on", &index, &[pair.clone()], true, None, false);
    assert_eq!(on.records.len(), 2);
    for (mate, (rec, row)) in on.records.iter().zip(&on.rows).enumerate() {
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
