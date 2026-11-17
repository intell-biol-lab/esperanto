//! Synthetic candidate-discovery regression; no human read or genotype fixture.
use esperanto_map::{
    align::{AlignConfig, Aligner},
    extend::CigarOp,
    fasta::parse_fasta_bytes,
    gtf::{Junction, JunctionLib},
    index::Index,
    seed::SeedParams,
};
use std::sync::Arc;

fn body(n: usize) -> Vec<u8> {
    let mut state = 741_u64;
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn fixture() -> (Index<'static>, Vec<u8>, JunctionLib) {
    let mut genome = body(1200);
    let first = genome[200..239].to_vec();
    let second = genome[551..563].to_vec();
    genome[239..241].copy_from_slice(b"GT");
    genome[549..551].copy_from_slice(b"AG");
    let core = [first, second].concat();
    let read = [vec![b'N'; 22], core.clone()].concat();
    let mut remote = body(200);
    remote[50..101].copy_from_slice(&core);
    remote[74] = if remote[74] == b'A' { b'C' } else { b'A' };
    let index = Index::build(
        parse_fasta_bytes(
            &[
                b">source\n".as_slice(),
                &genome,
                b"\n>remote\n",
                &remote,
                b"\n",
            ]
            .concat(),
        )
        .unwrap(),
        SeedParams::rna_default(),
    );
    let j = Junction {
        contig: 0,
        start: 239,
        end: 551,
        minus_strand: false,
    };
    (index, read, JunctionLib::build([j, j]))
}

#[test]
fn short_exon_library_candidate_competes_with_a_longer_remote_seed_chain() {
    let (index, read, lib) = fixture();
    let mut config = AlignConfig::rna_default();
    config.short_exon_candidates = true;
    let mut aligner = Aligner::new(&index, config).with_lib(Arc::new(lib));
    let observed = aligner.align_read(&read).unwrap();
    assert_eq!((observed.contig, observed.pos), (0, 200), "{observed:?}");
    assert_eq!(observed.second_score, 96, "remote competitor must remain");
    assert_eq!(
        observed.cigar,
        vec![
            CigarOp::SoftClip(22),
            CigarOp::Match(39),
            CigarOp::RefSkip(312),
            CigarOp::Match(12)
        ]
    );
}

#[test]
fn short_exon_supplement_is_default_off_and_requires_a_library() {
    let (index, read, lib) = fixture();
    let default = AlignConfig::rna_default();
    assert!(!default.short_exon_candidates);
    assert!(!AlignConfig::default().short_exon_candidates);
    let baseline = Aligner::new(&index, default)
        .with_lib(Arc::new(lib))
        .align_read(&read)
        .unwrap();
    assert_eq!((baseline.contig, baseline.pos), (1, 50));
    let mut config = default;
    config.short_exon_candidates = true;
    let no_library = Aligner::new(&index, config).align_read(&read).unwrap();
    assert_eq!(format!("{baseline:?}"), format!("{no_library:?}"));
    let unrelated = JunctionLib::build([Junction {
        contig: 0,
        start: 900,
        end: 1000,
        minus_strand: false,
    }]);
    let other_library = Aligner::new(&index, config)
        .with_lib(Arc::new(unrelated))
        .align_read(&read)
        .unwrap();
    assert_eq!(format!("{baseline:?}"), format!("{other_library:?}"));
}

#[test]
fn short_exon_supplement_preserves_reverse_query_coordinates() {
    let (index, mut read, lib) = fixture();
    read.reverse();
    for b in &mut read {
        *b = match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        };
    }
    let mut config = AlignConfig::rna_default();
    config.short_exon_candidates = true;
    let observed = Aligner::new(&index, config)
        .with_lib(Arc::new(lib))
        .align_read(&read)
        .unwrap();
    assert_eq!((observed.contig, observed.pos), (0, 200));
    assert_eq!(observed.strand, esperanto_map::seed::Strand::Minus);
    assert_eq!(
        observed.cigar,
        vec![
            CigarOp::SoftClip(22),
            CigarOp::Match(39),
            CigarOp::RefSkip(312),
            CigarOp::Match(12)
        ]
    );
}

#[test]
fn short_exon_supplement_does_not_replace_full_body_work() {
    let (index, _, lib) = fixture();
    let mut read = Vec::new();
    index.reference.contigs[0].decode_into(800, 875, &mut read);
    let lib = Arc::new(lib);
    let config = AlignConfig::rna_default();
    let before = Aligner::new(&index, config)
        .with_lib(lib.clone())
        .align_read(&read)
        .unwrap();
    let mut enabled = config;
    enabled.short_exon_candidates = true;
    let after = Aligner::new(&index, enabled)
        .with_lib(lib)
        .align_read(&read)
        .unwrap();
    assert_eq!(format!("{before:?}"), format!("{after:?}"));
}

#[test]
fn short_exon_supplement_preserves_left_tail_and_minus_transcript() {
    let (index, mut read, _) = fixture();
    let rc = |seq: &mut Vec<u8>| {
        seq.reverse();
        for b in seq {
            *b = match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            };
        }
    };
    let mut source = index.reference.contigs[0].slice_ascii(0, 1200);
    let mut remote = index.reference.contigs[1].slice_ascii(0, 200);
    rc(&mut source);
    rc(&mut remote);
    rc(&mut read);
    let mirrored = Index::build(
        parse_fasta_bytes(
            &[
                b">source\n".as_slice(),
                &source,
                b"\n>remote\n",
                &remote,
                b"\n",
            ]
            .concat(),
        )
        .unwrap(),
        SeedParams::rna_default(),
    );
    let j = Junction {
        contig: 0,
        start: 649,
        end: 961,
        minus_strand: true,
    };
    let lib = Arc::new(JunctionLib::build([j, j]));
    let mut config = AlignConfig::rna_default();
    config.short_exon_candidates = true;
    let observed = Aligner::new(&mirrored, config)
        .with_lib(lib)
        .align_read(&read)
        .unwrap();
    assert_eq!((observed.contig, observed.pos), (0, 637), "{observed:?}");
    assert_eq!(
        observed.cigar,
        vec![
            CigarOp::Match(12),
            CigarOp::RefSkip(312),
            CigarOp::Match(39),
            CigarOp::SoftClip(22)
        ]
    );
}

#[test]
fn short_exon_supplement_preserves_an_existing_high_coverage_genomic_path() {
    let (source, mut read, lib) = fixture();
    read[..22].copy_from_slice(&body(400)[350..372]);
    let mut remote = body(200);
    remote[50..123].copy_from_slice(&read);
    for offset in [42, 46, 50, 54, 58, 62, 66, 70] {
        let base = &mut remote[50 + offset];
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    let genome = source.reference.contigs[0].slice_ascii(0, 1200);
    let index = Index::build(
        parse_fasta_bytes(
            &[
                b">source\n".as_slice(),
                &genome,
                b"\n>remote\n",
                &remote,
                b"\n",
            ]
            .concat(),
        )
        .unwrap(),
        SeedParams::rna_default(),
    );
    let lib = Arc::new(lib);
    let config = AlignConfig::rna_default();
    let before = Aligner::new(&index, config)
        .with_lib(lib.clone())
        .align_read(&read)
        .unwrap();
    assert_eq!((before.contig, before.pos), (1, 50), "{before:?}");
    assert_eq!(
        before.cigar,
        vec![CigarOp::Match(70), CigarOp::SoftClip(3)],
        "{before:?}"
    );
    let mut enabled = config;
    enabled.short_exon_candidates = true;
    let after = Aligner::new(&index, enabled)
        .with_lib(lib)
        .align_read(&read)
        .unwrap();
    assert_eq!(format!("{before:?}"), format!("{after:?}"));
}

#[test]
fn short_exon_supplement_preserves_complete_l1_terminal_edit_on_both_strands() {
    let mut genome = body(1200);
    genome[239..241].copy_from_slice(b"GT");
    genome[549..551].copy_from_slice(b"AG");
    genome[561] = b'A';
    let mut read = [genome[200..239].to_vec(), genome[551..563].to_vec()].concat();
    read[49] = b'G';
    let directory =
        std::env::temp_dir().join(format!("esp-short-exon-l1-edit-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let bytes = [b">source\n".to_vec(), genome, b"\n".to_vec()].concat();
    let fasta = directory.join("reference.fa");
    std::fs::write(&fasta, &bytes).unwrap();
    std::fs::write(
        directory.join("reference.fa.fai"),
        b"source\t1200\t8\t1200\t1201\n",
    )
    .unwrap();
    let gtf = directory.join("annotation.gtf");
    std::fs::write(
        &gtf,
        concat!(
            "source\tx\ttranscript\t151\t700\t.\t+\t.\ttranscript_id \"a\";\n",
            "source\tx\texon\t151\t239\t.\t+\t.\ttranscript_id \"a\";\n",
            "source\tx\texon\t552\t700\t.\t+\t.\ttranscript_id \"a\";\n",
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
    let index = Index::build(
        parse_fasta_bytes(&bytes).unwrap(),
        SeedParams::rna_default(),
    );
    let lib = Arc::new(JunctionLib::build([Junction {
        contig: 0,
        start: 239,
        end: 551,
        minus_strand: false,
    }]));
    for reverse in [false, true] {
        let query: Vec<_> = if reverse {
            read.iter()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect()
        } else {
            read.clone()
        };
        let config = AlignConfig::rna_default();
        let before = Aligner::new(&index, config)
            .with_l1(l1.clone())
            .with_lib(lib.clone())
            .align_read(&query)
            .unwrap();
        assert_eq!(
            before.cigar,
            vec![
                CigarOp::Match(39),
                CigarOp::RefSkip(312),
                CigarOp::Match(12)
            ]
        );
        assert_eq!(before.ea_ag, 1);
        let mut enabled = config;
        enabled.short_exon_candidates = true;
        let after = Aligner::new(&index, enabled)
            .with_l1(l1.clone())
            .with_lib(lib.clone())
            .align_read(&query)
            .unwrap();
        assert_eq!(format!("{before:?}"), format!("{after:?}"));
    }
}
