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
