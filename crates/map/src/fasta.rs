//! Reference genome storage in 2-bit packed form.
//!
//! Each byte packs 4 bases, most-significant base first: base at `pos` lives
//! in byte `pos / 4` at shift `6 - 2 * (pos % 4)`. Tail bits of the final byte
//! (when `len % 4 != 0`) are zero. `N` cannot be encoded in 2 bits, so N runs
//! are kept as half-open intervals per contig and consulted on every decode.

use crate::error::AlignError;
use sha2::{Digest, Sha256};
use std::fs;
use std::path::Path;

/// A single nucleotide.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine.
    T,
    /// Any ambiguity code.
    N,
}

impl Base {
    /// Case-insensitive ASCII decoder; anything unrecognized becomes `N`.
    pub fn from_ascii(b: u8) -> Base {
        match b {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            _ => Base::N,
        }
    }

    /// 2-bit code; `None` for `N`.
    pub fn code(self) -> Option<u8> {
        match self {
            Base::A => Some(0),
            Base::C => Some(1),
            Base::G => Some(2),
            Base::T => Some(3),
            Base::N => None,
        }
    }
}

/// Half-open `[start, end)` interval of N bases.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Interval {
    /// First N position (0-based).
    pub start: u32,
    /// One past the last N position.
    pub end: u32,
}

/// One contig: name, length, 2-bit packed bases and N intervals.
pub struct Contig {
    /// Contig name (first whitespace-delimited token of the FASTA header).
    pub name: String,
    /// Sequence length in bases.
    pub len: u32,
    /// Packed bases; `ceil(len / 4)` bytes, high bases first within a byte.
    pub packed: Vec<u8>,
    /// Merged half-open N runs, sorted ascending.
    pub n_intervals: Vec<Interval>,
}

impl Contig {
    /// Decode `[start, end)` into ASCII bases, writing `N` over N intervals.
    ///
    /// Panics if the range is out of bounds; callers use validated ranges.
    pub fn decode_into(&self, start: u32, end: u32, out: &mut Vec<u8>) {
        let mut pos = start;
        while pos < end {
            let (is_n, stop) = self.ambiguity_span(pos, end);
            if is_n {
                let new_len = out
                    .len()
                    .checked_add((stop - pos) as usize)
                    .expect("capacity overflow");
                out.resize(new_len, b'N');
            } else {
                for p in pos..stop {
                    let byte = self.packed[(p >> 2) as usize];
                    let shift = 6 - 2 * (p & 3);
                    out.push(b"ACGT"[((byte >> shift) & 3) as usize]);
                }
            }
            pos = stop;
        }
    }

    // Return the N decision and the exclusive upper bound over which every
    // comparison on this exact is_n decision-tree path stays unchanged.
    // This also preserves old behavior for publicly constructible malformed
    // interval tables: do not assume sorting/merging beyond what is_n did.
    fn ambiguity_span(&self, pos: u32, end: u32) -> (bool, u32) {
        let (mut lo, mut hi, mut stop) = (0, self.n_intervals.len(), end);
        while lo < hi {
            #[cfg(test)]
            decode_tests::N_PROBES.with(|n| n.set(n.get() + 1));
            let mid = (lo + hi) / 2;
            let interval = &self.n_intervals[mid];
            if pos < interval.start {
                stop = stop.min(interval.start);
                hi = mid;
            } else if pos >= interval.end {
                lo = mid + 1;
            } else {
                return (true, stop.min(interval.end));
            }
        }
        (false, stop)
    }

    /// Append `[start, end)` ASCII bases to `out` (same as [`decode_into`]).
    pub fn decode_append(&self, start: u32, end: u32, out: &mut Vec<u8>) {
        self.decode_into(start, end, out);
    }

    /// ASCII slice of `[start, end)`; `N` wherever an N interval covers.
    pub fn slice_ascii(&self, start: u32, end: u32) -> Vec<u8> {
        let mut out = Vec::with_capacity((end - start) as usize);
        self.decode_into(start, end, &mut out);
        out
    }

    /// Bounds-checked base access: None when `pos` is past the contig end.
    pub fn base_checked(&self, pos: u32) -> Option<Base> {
        if pos >= self.len {
            None
        } else {
            Some(self.base(pos))
        }
    }

    /// Base at `pos`: N-interval table first, then 2-bit unpack.
    pub fn base(&self, pos: u32) -> Base {
        if self.is_n(pos) {
            return Base::N;
        }
        let byte = self.packed[(pos >> 2) as usize];
        let shift = 6 - 2 * (pos & 3);
        match (byte >> shift) & 3 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        }
    }

    fn is_n(&self, pos: u32) -> bool {
        let mut lo = 0usize;
        let mut hi = self.n_intervals.len();
        while lo < hi {
            #[cfg(test)]
            decode_tests::N_PROBES.with(|n| n.set(n.get() + 1));
            let mid = (lo + hi) / 2;
            let iv = &self.n_intervals[mid];
            if pos < iv.start {
                hi = mid;
            } else if pos >= iv.end {
                lo = mid + 1;
            } else {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod decode_tests {
    use super::{Contig, Interval};
    use std::cell::Cell;

    thread_local! {
        pub(super) static N_PROBES: Cell<usize> = const { Cell::new(0) };
    }

    fn contig(len: u32, n_intervals: Vec<Interval>) -> Contig {
        Contig {
            name: "decode-fixture".into(),
            len,
            packed: (0..len.div_ceil(4))
                .map(|i| (i.wrapping_mul(73).wrapping_add(27)) as u8)
                .collect(),
            n_intervals,
        }
    }

    // Independent scalar oracle; retains the original decision tree even for
    // arbitrary publicly constructible interval tables, not just valid FASTA.
    fn scalar(c: &Contig, start: u32, end: u32, out: &mut Vec<u8>) {
        for pos in start..end {
            let (mut lo, mut hi, mut is_n) = (0, c.n_intervals.len(), false);
            while lo < hi {
                let mid = (lo + hi) / 2;
                let interval = &c.n_intervals[mid];
                if pos < interval.start {
                    hi = mid;
                } else if pos >= interval.end {
                    lo = mid + 1;
                } else {
                    is_n = true;
                    break;
                }
            }
            out.push(if is_n {
                b'N'
            } else {
                b"ACGT"[((c.packed[(pos >> 2) as usize] >> (6 - 2 * (pos & 3))) & 3) as usize]
            });
        }
    }

    #[test]
    fn decode_ranges_preserve_prefix_and_scalar_ambiguity_semantics() {
        let tables = [
            vec![],
            vec![Interval { start: 0, end: 37 }],
            vec![
                Interval { start: 1, end: 3 },
                Interval { start: 3, end: 6 },
                Interval { start: 9, end: 11 },
                Interval { start: 35, end: 37 },
            ],
            // Public fields and old loader permit these: avoid silently
            // assuming merged or sorted interval tables in the optimization.
            vec![
                Interval { start: 17, end: 28 },
                Interval { start: 2, end: 5 },
                Interval { start: 8, end: 32 },
                Interval { start: 19, end: 12 },
            ],
        ];
        for intervals in tables {
            let c = contig(37, intervals);
            for start in 0..=37 {
                for end in start..=37 {
                    let mut expected = b"existing-prefix".to_vec();
                    scalar(&c, start, end, &mut expected);
                    let mut actual = b"existing-prefix".to_vec();
                    c.decode_into(start, end, &mut actual);
                    assert_eq!(actual, expected, "{start}..{end}");
                    let mut appended = b"existing-prefix".to_vec();
                    c.decode_append(start, end, &mut appended);
                    assert_eq!(appended, expected);
                }
            }
            let mut unchanged = b"prefix".to_vec();
            c.decode_into(20, 10, &mut unchanged);
            c.decode_into(u32::MAX, u32::MAX, &mut unchanged);
            assert_eq!(unchanged, b"prefix");
        }
    }

    #[test]
    fn decode_matches_scalar_for_arbitrary_interval_decision_trees() {
        let mut state = 0x68a93e7103d95b21u64;
        for _case in 0..256 {
            let mut intervals = Vec::new();
            for _ in 0..16 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let start = ((state >> 32) % 96) as u32;
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let end = ((state >> 32) % 96) as u32;
                intervals.push(Interval { start, end });
            }
            let c = contig(96, intervals);
            for start in 0..96 {
                let end = (start + 31).min(96);
                let mut expected = Vec::new();
                scalar(&c, start, end, &mut expected);
                assert_eq!(c.slice_ascii(start, end), expected);
            }
        }
    }

    #[test]
    fn decode_reuses_ambiguity_decisions_across_a_long_plain_span() {
        let c = contig(
            8192,
            vec![
                Interval { start: 0, end: 8 },
                Interval {
                    start: 8000,
                    end: 8010,
                },
            ],
        );
        let mut expected = Vec::new();
        scalar(&c, 64, 4160, &mut expected);
        N_PROBES.with(|n| n.set(0));
        let actual = c.slice_ascii(64, 4160);
        let probes = N_PROBES.with(Cell::get);
        assert_eq!(actual, expected);
        assert!(
            probes <= 4,
            "Repeated N-table probes for one plain span: {probes}"
        );
    }

    #[test]
    fn decode_preserves_ordinary_bounds_panics_and_ambiguity_first_access() {
        use std::panic::{catch_unwind, AssertUnwindSafe};
        for intervals in [
            vec![],
            vec![Interval { start: 2, end: 12 }],
            vec![Interval { start: 2, end: 7 }],
        ] {
            let mut c = contig(4, intervals);
            c.packed.truncate(1);
            for end in [4, 7, 12, 15] {
                let mut expected = b"prefix".to_vec();
                let mut actual = expected.clone();
                let old = catch_unwind(AssertUnwindSafe(|| scalar(&c, 0, end, &mut expected)));
                let new = catch_unwind(AssertUnwindSafe(|| c.decode_into(0, end, &mut actual)));
                assert_eq!(new.is_err(), old.is_err());
                assert_eq!(actual, expected);
            }
        }
    }
}

impl Base {
    pub fn to_ascii(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'N',
        }
    }
    /// Complement (A↔T, C↔G, N→N).
    pub fn complement(self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::N => Base::N,
        }
    }
}

/// A parsed reference: contigs plus the FASTA file's sha256.
pub struct Reference {
    /// Contigs in file order.
    pub contigs: Vec<Contig>,
    /// sha256 of the raw FASTA file bytes.
    pub fasta_sha256: [u8; 32],
}

impl Reference {
    /// Look up a contig by name; `None` if absent.
    pub fn contig_index(&self, name: &[u8]) -> Option<u32> {
        self.contigs
            .iter()
            .position(|c| c.name.as_bytes() == name)
            .map(|i| i as u32)
    }
}

/// Parse a FASTA file from disk.
///
/// Errors (with line numbers) on: empty file, empty or duplicate contig names.
/// Whitespace inside sequence lines is stripped; illegal characters become `N`.
pub fn parse_fasta(path: &Path) -> Result<Reference, AlignError> {
    let bytes = fs::read(path).map_err(|source| AlignError::FastaIo {
        path: path.display().to_string(),
        source,
    })?;
    parse_fasta_bytes(&bytes)
}

/// Parse a FASTA file from in-memory bytes; see [`parse_fasta`].
pub fn parse_fasta_bytes(bytes: &[u8]) -> Result<Reference, AlignError> {
    let mut contigs: Vec<Contig> = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_len: u32 = 0;
    let mut cur_packed: Vec<u8> = Vec::new();
    let mut cur_n: Vec<Interval> = Vec::new();
    let mut saw_any = false;

    let finish = |contigs: &mut Vec<Contig>,
                  name: Option<String>,
                  len: u32,
                  mut packed: Vec<u8>,
                  n_intervals: Vec<Interval>| {
        if let Some(name) = name {
            // pad the tail byte even when the final bases are all N
            packed.resize(len.div_ceil(4) as usize, 0);
            contigs.push(Contig {
                name,
                len,
                packed,
                n_intervals,
            });
        }
    };

    for (idx, raw_line) in split_lines(bytes).into_iter().enumerate() {
        let lineno = idx + 1;
        let line: Vec<u8> = raw_line
            .iter()
            .copied()
            .filter(|b| !b.is_ascii_whitespace())
            .collect();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            saw_any = true;
            finish(
                &mut contigs,
                cur_name.take(),
                cur_len,
                std::mem::take(&mut cur_packed),
                std::mem::take(&mut cur_n),
            );
            cur_len = 0;
            let header = &raw_line[1..];
            let token: Vec<u8> = header
                .iter()
                .copied()
                .take_while(|b| !b.is_ascii_whitespace())
                .collect();
            if token.is_empty() {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: "empty contig name".to_string(),
                });
            }
            cur_name = Some(String::from_utf8_lossy(&token).into_owned());
        } else {
            if cur_name.is_none() {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: "sequence before any '>' header".to_string(),
                });
            }
            for &b in &line {
                push_base(b, cur_len, &mut cur_packed, &mut cur_n);
                cur_len += 1;
            }
        }
    }

    if !saw_any {
        return Err(AlignError::FastaFormat {
            line: 0,
            msg: "empty fasta file".to_string(),
        });
    }
    finish(&mut contigs, cur_name, cur_len, cur_packed, cur_n);

    for (i, c) in contigs.iter().enumerate() {
        for d in contigs.iter().skip(i + 1) {
            if c.name == d.name {
                return Err(AlignError::FastaFormat {
                    line: 0,
                    msg: format!("duplicate contig name '{}'", c.name),
                });
            }
        }
    }

    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let fasta_sha256: [u8; 32] = hasher.finalize().into();

    Ok(Reference {
        contigs,
        fasta_sha256,
    })
}

fn push_base(b: u8, pos: u32, packed: &mut Vec<u8>, n_intervals: &mut Vec<Interval>) {
    let base = Base::from_ascii(b);
    match base.code() {
        Some(code) => {
            let byte_idx = (pos >> 2) as usize;
            while byte_idx >= packed.len() {
                packed.push(0);
            }
            let shift = 6 - 2 * (pos & 3);
            packed[byte_idx] |= code << shift;
        }
        None => match n_intervals.last_mut() {
            Some(iv) if iv.end == pos => iv.end += 1,
            _ => n_intervals.push(Interval {
                start: pos,
                end: pos + 1,
            }),
        },
    }
}

fn split_lines(bytes: &[u8]) -> Vec<&[u8]> {
    let mut lines = Vec::new();
    let mut start = 0usize;
    for (i, b) in bytes.iter().enumerate() {
        if *b == b'\n' {
            lines.push(&bytes[start..i]);
            start = i + 1;
        }
    }
    if start < bytes.len() {
        lines.push(&bytes[start..]);
    }
    lines
}
