//! K-mer encoding and the N-skip reseed stream.
//!
//! 2-bit codes: `A=0 C=1 G=2 T=3`, case-insensitive, with the 5'
//! base in the highest bits of the window. The canonical form of a k-mer
//! is `min(forward, revcomp)`. Any non-ACGT byte (N / IUPAC / other)
//! invalidates every window overlapping it; after an invalid byte the
//! stream reseeds from the byte after it (no phantom codes).

/// 2-bit code of one base byte, case-insensitive; `None` for any non-ACGT.
#[inline]
pub fn base_code(b: u8) -> Option<u64> {
    #[cfg(test)]
    tests::BASE_CODE_CALLS.with(|calls| calls.set(calls.get() + 1));
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a packed 2-bit k-mer code of length `k`
/// (complement every code and reverse their order).
pub fn revcomp_code(code: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for i in 0..k {
        let c = (code >> (2 * i)) & 3;
        rc = (rc << 2) | (3 - c);
    }
    rc
}

/// Canonical form of a forward k-mer code: `min(forward, revcomp)`.
pub fn canonical(code: u64, k: usize) -> u64 {
    code.min(revcomp_code(code, k))
}

/// Forward 2-bit code of `seq[off..off + k]`; `None` when the window is
/// out of bounds or contains a non-ACGT byte.
pub fn window_code(seq: &[u8], off: usize, k: usize) -> Option<u64> {
    if off + k > seq.len() {
        return None;
    }
    let mut code = 0u64;
    for &b in &seq[off..off + k] {
        code = (code << 2) | base_code(b)?;
    }
    Some(code)
}

/// Case-preserving reverse complement: `a<->t`, `c<->g`; non-ACGT bytes
/// pass through unchanged. Used for minus-strand oriented reads — k-mer
/// encoding is case-insensitive, this is consistency only.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Streaming forward k-mer codes over a read with N-skip reseed.
///
/// Yields `(offset, forward_code)` for every ACGT-only window in ascending
/// offset order. When a window contains a non-ACGT byte, the stream
/// reseeds at the byte after the first invalid byte in that window.
/// Each input byte is decoded at most once, retaining a rolling forward code.
#[derive(Debug, Clone)]
pub struct KmerStream<'a> {
    seq: &'a [u8],
    k: usize,
    // Cursor of the next input byte, rather than the next window start.
    pos: usize,
    code: u64,
    mask: u64,
    valid_bases: usize,
}

impl<'a> KmerStream<'a> {
    /// Create the stream; `k` must be in `1..=32`.
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!((1..=32).contains(&k), "k-mer size must be in 1..=32");
        KmerStream {
            seq,
            k,
            pos: 0,
            code: 0,
            // k=32 keeps all 64 bits without ever shifting by 64.
            mask: u64::MAX >> (2 * (32 - k)),
            valid_bases: 0,
        }
    }
}

impl<'a> Iterator for KmerStream<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<(usize, u64)> {
        while self.pos < self.seq.len() {
            let byte = self.seq[self.pos];
            self.pos += 1;
            match base_code(byte) {
                Some(bits) => {
                    self.code = ((self.code << 2) | bits) & self.mask;
                    if self.valid_bases < self.k {
                        self.valid_bases += 1;
                    }
                    if self.valid_bases == self.k {
                        return Some((self.pos - self.k, self.code));
                    }
                }
                None => {
                    self.code = 0;
                    self.valid_bases = 0;
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::KmerStream;
    use std::cell::Cell;

    thread_local! {
        // Count real decoding work on this test's thread; absent in production.
        pub(super) static BASE_CODE_CALLS: Cell<usize> = const { Cell::new(0) };
    }

    // Deliberately enumerate every window and encode it independently of
    // base_code/window_code and the stream's invalid-byte reseed strategy.
    fn slow_windows(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
        seq.windows(k)
            .enumerate()
            .filter_map(|(offset, window)| {
                let mut code = 0u64;
                for &byte in window {
                    let bits = match byte {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => return None,
                    };
                    code = (code << 2) | bits;
                }
                Some((offset, code))
            })
            .collect()
    }

    fn assert_matches_slow(seq: &[u8], k: usize) {
        assert_eq!(
            KmerStream::new(seq, k).collect::<Vec<_>>(),
            slow_windows(seq, k),
            "length={}, k={k}, prefix={:?}",
            seq.len(),
            &seq[..seq.len().min(80)]
        );
    }

    #[test]
    fn stream_matches_slow_exhaustive_tiny_alphabet() {
        const ALPHABET: &[u8] = b"AcGtN\xff";
        for len in 0..=6u32 {
            let mut seq = vec![0; len as usize];
            for mut word in 0..ALPHABET.len().pow(len) {
                for byte in &mut seq {
                    *byte = ALPHABET[word % ALPHABET.len()];
                    word /= ALPHABET.len();
                }
                for k in 1..=32 {
                    assert_matches_slow(&seq, k);
                }
            }
        }
    }

    #[test]
    fn stream_matches_slow_every_byte_at_window_boundaries() {
        for byte in 0..=u8::MAX {
            let mut seq = b"AcGt".repeat(10);
            seq.push(byte);
            seq.extend_from_slice(&b"tGcA".repeat(10));
            for k in 1..=32 {
                assert_matches_slow(&seq, k);
                assert_matches_slow(&seq[40..], k);
                assert_matches_slow(&seq[..41], k);
            }
        }
    }

    #[test]
    fn stream_matches_slow_long_deterministic_random_all_k() {
        // Fixed, local PRNG: exercise raw bytes as well as long mixed-case
        // valid runs, so large k values produce substantial output too.
        let mut state = 0x6a09_e667_f3bc_c909u64;
        let mut raw: Vec<u8> = (0..=u8::MAX).collect();
        let mut mixed = Vec::new();
        for i in 0..65_536 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let byte = (state >> 32) as u8;
            raw.push(byte);
            mixed.push(if i % 521 < 384 {
                b"ACGTacgt"[usize::from(byte) % 8]
            } else {
                byte
            });
        }
        for k in 1..=32 {
            assert_matches_slow(&raw, k);
            assert_matches_slow(&mixed, k);
        }
    }

    #[test]
    fn stream_clone_resumes_independently_including_exhaustion() {
        let mut seq = b"AcGt".repeat(12);
        seq.extend_from_slice(b"NNry\x00\xff");
        seq.extend_from_slice(&b"tGcA".repeat(12));
        seq.extend_from_slice(b"Nacg");
        for k in 1..=32 {
            let expected = slow_windows(&seq, k);
            let mut stream = KmerStream::new(&seq, k);
            for consumed in 0..=expected.len() {
                assert_eq!(
                    stream.clone().collect::<Vec<_>>(),
                    expected[consumed..],
                    "k={k}, consumed={consumed}"
                );
                assert_eq!(stream.next(), expected.get(consumed).copied());
            }
            for _ in 0..3 {
                assert_eq!(stream.next(), None);
                assert_eq!(stream.clone().next(), None);
            }
            let mut short = KmerStream::new(&seq[..k - 1], k);
            assert_eq!(short.clone().next(), None);
            assert_eq!(short.next(), None);
            assert_eq!(short.clone().next(), None);
        }
    }

    #[test]
    fn stream_k32_keeps_all_bits_and_reseeds_without_overflow() {
        let mut seq = vec![b'T'; 64];
        seq.push(b'N');
        seq.extend_from_slice(&b"acgt".repeat(16));
        assert_matches_slow(&seq, 32);
        assert_eq!(KmerStream::new(&seq, 32).next(), Some((0, u64::MAX)));
        for invalid_at in 0..96 {
            let mut interrupted = vec![b't'; 96];
            interrupted[invalid_at] = b'R';
            assert_matches_slow(&interrupted, 32);
        }
        for len in 0..32 {
            assert_matches_slow(&seq[..len], 32);
        }
    }

    #[test]
    fn stream_decodes_each_byte_at_most_once() {
        let clean = b"AcGt".repeat(1024);
        let mut interrupted = clean.clone();
        for i in (47..interrupted.len()).step_by(97) {
            interrupted[i] = [b'N', b'r', 0, 255][i % 4];
        }
        let mut clean_counts = Vec::new();
        let mut interrupted_counts = Vec::new();
        for (seq, counts) in [
            (&clean, &mut clean_counts),
            (&interrupted, &mut interrupted_counts),
        ] {
            for k in 1..=32 {
                BASE_CODE_CALLS.with(|calls| calls.set(0));
                let actual = KmerStream::new(seq, k).collect::<Vec<_>>();
                counts.push(BASE_CODE_CALLS.with(Cell::get));
                assert_eq!(actual, slow_windows(seq, k));
            }
        }
        eprintln!(
            "base_code calls for {} bytes, k=1..32: clean={clean_counts:?}; interrupted={interrupted_counts:?}",
            clean.len()
        );
        assert_eq!(
            clean_counts,
            vec![clean.len(); 32],
            "overlapping valid windows must reuse the rolling code"
        );
        assert!(
            interrupted_counts
                .iter()
                .all(|&calls| calls <= interrupted.len()),
            "invalid-byte reseeding must not decode a byte twice"
        );
    }
}
