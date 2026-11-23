use super::*;
use std::cell::Cell;
use std::collections::hash_map::{DefaultHasher, RandomState};
use std::hash::{BuildHasher, Hasher};
use std::rc::Rc;

struct CountBuild {
    state: RandomState,
    calls: Rc<Cell<usize>>,
}

struct CountHash {
    inner: DefaultHasher,
    calls: Rc<Cell<usize>>,
}

impl Hasher for CountHash {
    fn finish(&self) -> u64 {
        self.calls.set(self.calls.get() + 1);
        self.inner.finish()
    }

    fn write(&mut self, bytes: &[u8]) {
        self.inner.write(bytes);
    }
}

impl BuildHasher for CountBuild {
    type Hasher = CountHash;

    fn build_hasher(&self) -> Self::Hasher {
        CountHash {
            inner: self.state.build_hasher(),
            calls: self.calls.clone(),
        }
    }
}

fn observation(base: u8, bq: u8) -> Observation {
    Observation {
        base,
        bq,
        mapq: 60.0,
        end_fraction: 0.5,
        orientation: 0,
        softclip_fraction: 0.0,
        gap_nearness: 0.0,
    }
}

#[test]
fn template_name_borrows_exact_normalized_bytes() {
    for input in [
        &b""[..],
        b"/1",
        b"/2",
        b"a/1",
        b"a/2",
        b"a/3",
        b"a/12",
        b"a/1x",
        b"a/1/2",
        b"a",
        b"\xff/1",
    ] {
        let actual = template_name(input);
        let expected = if input.ends_with(b"/1") || input.ends_with(b"/2") {
            &input[..input.len() - 2]
        } else {
            input
        };
        assert_eq!(&actual[..], expected);
        assert_eq!(
            actual.as_ptr(),
            input.as_ptr(),
            "normalization must borrow the original name without allocation"
        );
    }
}

#[test]
fn template_merge_hashes_once_and_preserves_conflict_recovery() {
    let calls = Rc::new(Cell::new(0));
    let mut entries = HashMap::with_capacity_and_hasher(
        16,
        CountBuild {
            state: RandomState::new(),
            calls: calls.clone(),
        },
    );
    entries.insert(1, observation(b'A', 30));
    for (base, bq, expected) in [
        (b'A', 30, b'A'),
        (b'G', 30, b'N'),
        (b'T', 20, b'N'),
        (b'G', 40, b'G'),
    ] {
        calls.set(0);
        merge_template(&mut entries, 2, observation(base, bq));
        assert_eq!(calls.get(), 1, "each insert/merge must hash its key once");
        assert_eq!(entries[&2].base, expected);
    }
    assert_eq!(entries[&2].bq, 40);
    assert_eq!(entries[&1].base, b'A');
    assert_eq!(entries.len(), 2);
}
