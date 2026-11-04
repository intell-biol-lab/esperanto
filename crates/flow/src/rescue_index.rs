//! One immutable transformed index, reused only across adjacent passes.

use esperanto_map::{error::AlignError, index::Index, index_io::OwnedIndex};
use std::path::{Path, PathBuf};

#[derive(Default)]
pub(crate) struct RescueIndex {
    current: Option<(PathBuf, OwnedIndex)>,
}

impl RescueIndex {
    pub(crate) fn index_for(&mut self, path: &Path) -> Result<Index<'_>, AlignError> {
        if self.current.as_ref().map(|(p, _)| p.as_path()) != Some(path) {
            // Drop the previous mapping before opening/decoding the next one,
            // including when that next open fails.
            self.current = None;
            self.current = Some((path.to_owned(), OwnedIndex::open(path)?));
        }
        Ok(self.current.as_ref().unwrap().1.index())
    }

    /// Release before verification when no following pass needs this index.
    pub(crate) fn release_unless(&mut self, next: Option<&Path>) {
        if self.current.as_ref().map(|(path, _)| path.as_path()) != next {
            self.current = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use esperanto_map::{fasta, index_io, seed::SeedParams};
    use std::{
        fs,
        sync::atomic::{AtomicU64, Ordering},
    };

    static NEXT: AtomicU64 = AtomicU64::new(0);

    struct Fixture(PathBuf);

    impl Fixture {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "esp-adjacent-rescue-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn index(&self, name: &str) -> PathBuf {
            let path = self.0.join(name);
            let reference = fasta::parse_fasta_bytes(b">chr1\nACGTACGTTGCAAGCTACGTACGT\n").unwrap();
            let index = Index::build(reference, SeedParams { k: 5, w: 2 });
            index_io::save(&index, &path).unwrap();
            path
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    fn adjacent_pass_keeps_existing_backing_even_if_path_is_unavailable() {
        let f = Fixture::new();
        let tc = f.index("tc.paidx");
        let mut indexes = RescueIndex::default();
        let (key, positions) = {
            let idx = indexes.index_for(&tc).unwrap();
            (idx.kmers[0], idx.query(idx.kmers[0]).positions.to_vec())
        };
        indexes.release_unless(Some(&tc));
        // A repeated open would fail here. The immutable mmap is still valid.
        fs::rename(&tc, f.0.join("hidden.paidx")).unwrap();
        let idx = indexes.index_for(&tc).expect("adjacent TC must not reopen");
        assert_eq!(idx.query(key).positions, positions);
        indexes.release_unless(None);
        assert!(indexes.current.is_none());
    }

    #[test]
    fn changed_path_releases_owner_before_open_failure() {
        let f = Fixture::new();
        let ag = f.index("ag.paidx");
        let mut indexes = RescueIndex::default();
        indexes.index_for(&ag).unwrap();
        assert!(indexes.index_for(&f.0.join("missing.paidx")).is_err());
        assert!(indexes.current.is_none());
        #[cfg(target_os = "linux")]
        assert!(!fs::read_to_string("/proc/self/maps")
            .unwrap()
            .contains(ag.to_str().unwrap()));
    }

    #[test]
    fn leaving_adjacent_run_unmaps_before_next_index_and_on_drop() {
        let f = Fixture::new();
        let ag = f.index("ag.paidx");
        let tc = f.index("tc.paidx");
        {
            let mut indexes = RescueIndex::default();
            indexes.index_for(&ag).unwrap();
            indexes.release_unless(Some(&tc));
            assert!(indexes.current.is_none());
            #[cfg(target_os = "linux")]
            assert!(!fs::read_to_string("/proc/self/maps")
                .unwrap()
                .contains(ag.to_str().unwrap()));
            indexes.index_for(&tc).unwrap();
            indexes.release_unless(Some(&tc));
            assert!(indexes.current.is_some());
        }
        #[cfg(target_os = "linux")]
        assert!(!fs::read_to_string("/proc/self/maps")
            .unwrap()
            .contains(tc.to_str().unwrap()));
    }
}
