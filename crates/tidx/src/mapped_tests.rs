use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);

fn fixture(version: u32, aligned: bool) -> Vec<u8> {
    let k = if version == 1 { 31 } else { 15 };
    let name = if aligned {
        b"tx0_____".as_slice()
    } else {
        b"tx0".as_slice()
    };
    let mut out = vec![0u8; HEADER_LEN];
    for (off, value) in [(4, version), (8, k), (12, M_FIELD), (16, 2), (64, 1)] {
        out[off..off + 4].copy_from_slice(&value.to_le_bytes());
    }
    out[..4].copy_from_slice(&MAGIC);
    for (off, value) in [
        (24, SEED),
        (80, 5),
        (88, 4),
        (96, name.len() as u64),
        (104, 2),
    ] {
        out[off..off + 8].copy_from_slice(&value.to_le_bytes());
    }
    out.extend_from_slice(name);
    for offset in [0u32, name.len() as u32] {
        out.extend_from_slice(&offset.to_le_bytes());
    }
    for offset in [0u64, 3, 3, 4, 5] {
        out.extend_from_slice(&offset.to_le_bytes());
    }
    for key in [5u64, 5, 9, 6, 7] {
        if version == 1 {
            out.extend_from_slice(&key.to_le_bytes());
        } else {
            out.extend_from_slice(&(key as u32).to_le_bytes());
        }
    }
    for offset in [1u32, 2, 3, 4, 0xfedcba98] {
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&offset.to_le_bytes());
    }
    out
}

fn temporary(bytes: &[u8]) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "esperanto-mapped-tidx-{}-{}.tidx",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::write(&path, bytes).unwrap();
    path
}

#[test]
fn mapped_hits_match_owning_reader_both_versions_and_unaligned_tables() {
    for version in [1, 2] {
        for aligned in [false, true] {
            let bytes = fixture(version, aligned);
            let old = parse_index(&bytes).unwrap();
            let path = temporary(&bytes);
            let mapped = MappedTidx::open(&path).unwrap();
            assert_eq!(mapped.k(), old.k());
            assert_eq!(mapped.tx_count(), old.tx_count());
            assert_eq!(mapped.transcript_name(0), old.transcript_name(0));
            assert_eq!(mapped.transcript_name(1), "");
            let shift = 2 * mapped.k() - 2;
            for bucket in 0..5u64 {
                for low in 0..12u64 {
                    let key = (bucket << shift) | low;
                    let hits = mapped.lookup_hits(key);
                    assert_eq!(hits.len(), old.lookup(key).len());
                    assert_eq!(hits.is_empty(), old.lookup(key).is_empty());
                    assert_eq!(hits.iter().collect::<Vec<_>>(), old.lookup(key));
                }
            }
            // Normal mapped access must not instantiate the legacy payload vector.
            assert!(mapped.legacy_payloads.get().is_none());
            assert_eq!(mapped.lookup(5), &[(0, 1), (0, 2)]);
            assert!(mapped.legacy_payloads.get().is_some());
            fs::remove_file(path).unwrap();
        }
    }
}

#[test]
fn mapped_v1_suffixes_preserve_bits_above_u32() {
    let high = (1u64 << 40) | 5;
    for aligned in [false, true] {
        let mut bytes = fixture(1, aligned);
        let key_start = HEADER_LEN + if aligned { 8 } else { 3 } + 8 + 40;
        for (i, key) in [high, high, high + 4, high + 1, high + 2]
            .into_iter()
            .enumerate()
        {
            bytes[key_start + i * 8..key_start + (i + 1) * 8].copy_from_slice(&key.to_le_bytes());
        }
        let old = parse_index(&bytes).unwrap();
        let path = temporary(&bytes);
        let mapped = MappedTidx::open(&path).unwrap();
        for bucket in 0..4u64 {
            for low in [
                5,
                high - 1,
                high,
                high + 1,
                high + 2,
                high + 3,
                high + 4,
                high + 5,
            ] {
                let key = (bucket << 60) | low;
                assert_eq!(
                    mapped.lookup_hits(key).iter().collect::<Vec<_>>(),
                    old.lookup(key)
                );
            }
        }
        assert_eq!(
            mapped.lookup_hits(high).iter().collect::<Vec<_>>(),
            [(0, 1), (0, 2)]
        );
        assert!(mapped.legacy_payloads.get().is_none());
        fs::remove_file(path).unwrap();
    }
}

#[test]
fn mapped_clone_owns_mapping_after_original_and_file_handle_are_dropped() {
    let path = temporary(&fixture(1, false));
    let original = MappedTidx::open(&path).unwrap();
    let clone = original.clone();
    drop(original);
    // Unlinking does not mutate the mapped inode; the clone still owns its pages.
    fs::remove_file(path).unwrap();
    assert_eq!(
        clone.lookup_hits(5).iter().collect::<Vec<_>>(),
        [(0, 1), (0, 2)]
    );
}

#[test]
fn mapped_queries_are_shared_without_legacy_materialization() {
    let path = temporary(&fixture(2, false));
    let mapped = MappedTidx::open(&path).unwrap();
    let clone = mapped.clone();
    std::thread::scope(|scope| {
        for reader in [&mapped, &clone] {
            scope.spawn(move || {
                for _ in 0..1000 {
                    assert_eq!(
                        reader.lookup_hits(5).iter().collect::<Vec<_>>(),
                        [(0, 1), (0, 2)]
                    );
                }
            });
        }
    });
    assert!(mapped.legacy_payloads.get().is_none());
    assert!(clone.legacy_payloads.get().is_none());
    assert_eq!(clone.lookup(5), [(0, 1), (0, 2)]);
    assert!(mapped.legacy_payloads.get().is_some());
    fs::remove_file(path).unwrap();
}

#[test]
fn mapped_full_width_keys_and_legacy_view_iterator_preserve_values() {
    let mut bytes = fixture(1, true);
    bytes[8..12].copy_from_slice(&32u32.to_le_bytes());
    let path = temporary(&bytes);
    let old = parse_index(&bytes).unwrap();
    let mapped = MappedTidx::open(&path).unwrap();
    for key in [0, 5, 9, (2u64 << 62) | 6, (3u64 << 62) | 7, u64::MAX] {
        assert_eq!(
            mapped.lookup_hits(key).iter().collect::<Vec<_>>(),
            old.lookup(key)
        );
        let hits = LookupHits::from_slice(old.lookup(key));
        assert_eq!(hits.iter().len(), hits.len());
        assert_eq!(hits.iter().collect::<Vec<_>>(), old.lookup(key));
    }
    assert!(mapped.legacy_payloads.get().is_none());
    fs::remove_file(path).unwrap();
}

#[test]
fn mapped_reader_preserves_existing_format_rejections() {
    let clean = fixture(2, false);
    let mut cases = Vec::new();
    for end in [0, 4, 127, 128, clean.len() - 1] {
        cases.push(clean[..end].to_vec());
    }
    let mut trailing = clean.clone();
    trailing.push(0);
    cases.push(trailing);
    for (offset, value) in [
        (0, 0xff),
        (4, 99),
        (8, 0),
        (12, 0),
        (16, 0),
        (24, 0),
        (88, 0),
        (104, 0),
        (131, 1),
        (139, 1),
    ] {
        let mut corrupt = clean.clone();
        corrupt[offset] = value;
        cases.push(corrupt);
    }
    for offset in [80, 88, 96, 104] {
        let mut corrupt = clean.clone();
        corrupt[offset..offset + 8].fill(0xff);
        cases.push(corrupt);
    }
    for bytes in cases {
        let old_error = parse_index(&bytes).unwrap_err().to_string();
        let path = temporary(&bytes);
        let error = MappedTidx::open(&path).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), old_error);
        fs::remove_file(path).unwrap();
    }
}
