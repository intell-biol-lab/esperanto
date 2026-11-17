use esperanto_map::{fasta, index::Index, index_io, seed::SeedParams};
use std::{fs, path::PathBuf};

// Legacy process-lifetime storage allows hits to outlive the Index value.
fn legacy_static_hit(index: Index<'static>, key: u64) -> &'static [u64] {
    index.query(key).positions
}

#[test]
fn hits_borrow_backing_storage_instead_of_temporary_index_view() {
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGTACGTTGCAAGCT\n").unwrap();
    let built = Index::build(reference, SeedParams { k: 5, w: 2 });
    let key = built.kmers[0];
    let expected = built.query(key).positions.to_vec();
    assert_eq!(legacy_static_hit(built, key), expected);
}

#[test]
fn scoped_index_exact_hits_parallel_borrow_and_unmap() {
    let path = std::env::temp_dir().join(format!("esp-owned-paidx-{}.paidx", std::process::id()));
    let reference =
        fasta::parse_fasta_bytes(b">chr1\nACGTACGTTGCAAGCTACGTACGTACGTTGCAAGCT\n").unwrap();
    let built = Index::build(reference, SeedParams { k: 5, w: 2 });
    index_io::save(&built, &path).unwrap();
    {
        let owner = index_io::OwnedIndex::open(&path).unwrap();
        let view = owner.index();
        assert_eq!(view.kmers, built.kmers);
        assert_eq!(view.offsets, built.offsets);
        assert_eq!(view.counts, built.counts);
        assert_eq!(view.positions, built.positions);
        assert_eq!(view.reference.fasta_sha256, built.reference.fasta_sha256);
        std::thread::scope(|s| {
            for _ in 0..8 {
                s.spawn(|| {
                    for &key in built.kmers {
                        let a = built.query(key);
                        let b = view.query(key);
                        assert_eq!(a.count, b.count);
                        assert_eq!(a.positions, b.positions);
                    }
                    assert_eq!(view.query(u64::MAX).count, 0);
                });
            }
        });
        #[cfg(target_os = "linux")]
        assert!(mapped(&path));
    }
    #[cfg(target_os = "linux")]
    assert!(
        !mapped(&path),
        "directional index mmap remains after owner drops"
    );
    fs::remove_file(path).unwrap();
}

#[cfg(target_os = "linux")]
fn mapped(path: &PathBuf) -> bool {
    fs::read_to_string("/proc/self/maps")
        .unwrap()
        .contains(path.to_str().unwrap())
}

#[test]
fn scoped_reader_matches_legacy_rejections_and_cleans_failed_maps() {
    let path = std::env::temp_dir().join(format!("esp-bad-paidx-{}.paidx", std::process::id()));
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGTACGTTGCAAGCT\n").unwrap();
    let built = Index::build(reference, SeedParams { k: 5, w: 2 });
    let bytes = index_io::encode(&built).unwrap();
    let mut cases = vec![Vec::new(), bytes[..30].to_vec()];
    let mut bad = bytes.clone();
    bad[0] ^= 1;
    cases.push(bad);
    let mut bad = bytes.clone();
    bad[8..12].copy_from_slice(&99u32.to_le_bytes());
    cases.push(bad);
    let mut bad = bytes.clone();
    bad.push(1);
    cases.push(bad);
    let mut bad = bytes.clone();
    bad[28..36].copy_from_slice(&u64::MAX.to_le_bytes());
    cases.push(bad);
    let c = &built.reference.contigs[0];
    let contig_bytes = 4 + c.name.len() + 4 + 4 + c.packed.len() + 4;
    let pad = (8 - (60 + contig_bytes) % 8) % 8;
    let kmers_start = 76 + contig_bytes + pad;
    assert!(built.kmers.len() >= 2);
    let mut bad = bytes.clone();
    let first = bad[kmers_start..kmers_start + 8].to_vec();
    bad[kmers_start + 8..kmers_start + 16].copy_from_slice(&first);
    cases.push(bad);
    let mut bad = bytes.clone();
    let offsets_start = kmers_start + built.kmers.len() * 8;
    bad[offsets_start..offsets_start + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    cases.push(bad);
    for bytes in cases {
        fs::write(&path, &bytes).unwrap();
        assert!(index_io::OwnedIndex::open(&path).is_err());
        // Empty mmap returns an I/O error; compare legacy decoder for all
        // nonempty files, including exact version error paths.
        if !bytes.is_empty() {
            let old = index_io::decode(Box::leak(bytes.into_boxed_slice()), path.to_str().unwrap())
                .err()
                .unwrap()
                .to_string();
            let new = index_io::OwnedIndex::open(&path).err().unwrap().to_string();
            assert_eq!(old, new);
        }
        #[cfg(target_os = "linux")]
        assert!(!mapped(&path), "failed load leaked a mapping");
    }
    fs::remove_file(path).unwrap();
}
