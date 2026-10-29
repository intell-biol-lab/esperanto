//! Minimizer seeding (RNA high-density parameters).
//!
//! K-mers are encoded in 2-bit canonical form (min of forward and reverse
//! complement); window sort keys are the splitmix64 finalizer of the canonical
//! k-mer (not the raw k-mer value). K-mers overlapping an N (code 0xFF) are
//! invalid; a window with invalid members still selects among its valid ones.
//! Ties pick the leftmost k-mer; a minimizer shared by adjacent windows is
//! emitted only once.

use crate::fasta::Base;
use std::collections::VecDeque;

/// Minimizer extraction parameters.
#[derive(Clone, Copy, Debug)]
pub struct SeedParams {
    /// K-mer length in bases.
    pub k: u32,
    /// Window size in k-mers.
    pub w: u32,
}

impl SeedParams {
    /// RNA defaults: k=15, w=5.
    pub fn rna_default() -> SeedParams {
        SeedParams { k: 15, w: 5 }
    }
}

/// Strand of a seed relative to its source sequence.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Strand {
    /// Forward orientation.
    #[default]
    Plus,
    /// Reverse-complement orientation.
    Minus,
}

impl Strand {
    /// Opposite strand.
    pub fn flip(self) -> Strand {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }
}

/// One minimizer: canonical 2-bit k-mer at a position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Minimizer {
    /// Canonical (min of k-mer and revcomp) 2-bit encoding.
    pub kmer: u64,
    /// 0-based position of the k-mer start.
    pub pos: u32,
    /// Orientation of the original k-mer (Plus if it equals the canonical).
    pub strand: Strand,
}

/// A seed hit pairing a read minimizer with a reference position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Anchor {
    /// Read position (minimizer start).
    pub qpos: u32,
    /// Reference position (minimizer start).
    pub rpos: u32,
    /// Contig index.
    pub contig: u32,
    /// Seed strand (read-relative XOR reference-relative).
    pub strand: Strand,
}

/// splitmix64 finalizer — window sort key for a canonical k-mer.
#[inline]
pub fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

const N_CODE: u64 = 0xFF;

/** Reusable buffers for minimizer extraction (per-thread scratch).
 *
 * Allocation-elimination only: every buffer is cleared before refill, so
 * results are identical to the fresh-`Vec` path.
 */
#[derive(Default)]
pub struct SeedScratch {
    codes: Vec<u64>,
    keys: Vec<Option<(u64, u64, Strand)>>,
    dq: VecDeque<usize>,
}

/// Extract minimizers from an ASCII read sequence.
pub fn minimizers(seq: &[u8], params: SeedParams) -> Vec<Minimizer> {
    let mut out = Vec::new();
    minimizers_into(seq, params, &mut SeedScratch::default(), &mut out);
    out
}

/// [`minimizers`] into reusable buffers (`out` is cleared, then filled).
pub fn minimizers_into(
    seq: &[u8],
    params: SeedParams,
    scratch: &mut SeedScratch,
    out: &mut Vec<Minimizer>,
) {
    let SeedScratch { codes, keys, dq } = scratch;
    codes.clear();
    codes.extend(seq.iter().map(|&b| match Base::from_ascii(b).code() {
        Some(c) => c as u64,
        None => N_CODE,
    }));
    minimizers_core(codes, params.k, params.w, keys, dq, out);
}

/// Complete the query's two incomplete minimizer-window edges.
///
/// Call immediately after `minimizers_into` for the same query and parameters.
/// Only valid, canonical k-mers from the first/last `w - 1` starts are added;
/// reference sketching and the original minimizer extraction are unchanged.
/// The caller must apply the ordinary occurrence/anchor budgets to the union.
pub(crate) fn append_query_edge_kmers(
    read_len: usize,
    params: SeedParams,
    scratch: &SeedScratch,
    out: &mut Vec<Minimizer>,
) -> usize {
    let k = params.k as usize;
    if k == 0 || k > 32 || params.w <= 1 || read_len < k {
        return 0;
    }
    let n = read_len - k + 1;
    if scratch.codes.len() != read_len || scratch.keys.len() != n || n > u32::MAX as usize {
        return 0;
    }
    let before = out.len();
    let edge = (params.w as usize - 1).min(n);
    for pos in (0..edge).chain(n.saturating_sub(edge).max(edge)..n) {
        if let Some((_, kmer, strand)) = scratch.keys[pos] {
            out.push(Minimizer {
                kmer,
                pos: pos as u32,
                strand,
            });
        }
    }
    out.sort_unstable_by_key(|m| (m.pos, m.kmer));
    out.dedup_by_key(|m| (m.pos, m.kmer));
    out.len() - before
}

/** Extract minimizers from pre-encoded 2-bit codes (`N` must be 0xFF).
 *
 * Shared by read seeding and index construction. Uses a monotonic queue so
 * the whole pass is O(n).
 */
pub fn minimizers_from_codes(codes: &[u64], k: u32, w: u32) -> Vec<Minimizer> {
    let mut out = Vec::new();
    minimizers_core(codes, k, w, &mut Vec::new(), &mut VecDeque::new(), &mut out);
    out
}

fn minimizers_core(
    codes: &[u64],
    k: u32,
    w: u32,
    keys: &mut Vec<Option<(u64, u64, Strand)>>,
    dq: &mut VecDeque<usize>,
    out: &mut Vec<Minimizer>,
) {
    out.clear();
    let k = k as usize;
    let w = w as usize;
    if codes.len() < k || w == 0 || k == 0 || k > 32 {
        return;
    }
    let n_kmers = codes.len() - k + 1;
    // Per-k-mer (sort key, kmer, strand); k-mers overlapping an N are invalid.
    keys.clear();
    keys.reserve(n_kmers);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut valid_run = 0usize; // current N-free run length
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    for (i, &c) in codes.iter().enumerate() {
        if c == N_CODE {
            valid_run = 0;
            fwd = 0;
            rev = 0;
        } else {
            fwd = ((fwd << 2) | c) & mask;
            // revcomp: new base complemented into the high slot
            rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
            valid_run += 1;
        }
        if i + 1 >= k {
            if valid_run >= k {
                let (kmer, strand) = if fwd <= rev {
                    (fwd, Strand::Plus)
                } else {
                    (rev, Strand::Minus)
                };
                keys.push(Some((mix64(kmer), kmer, strand)));
            } else {
                keys.push(None);
            }
        }
    }
    debug_assert_eq!(keys.len(), n_kmers);

    // Sliding-window minimum (monotonic queue); a window containing invalid
    // k-mers still selects among its VALID members (only fully-invalid
    // windows emit nothing).
    dq.clear();
    for i in 0..n_kmers {
        if keys[i].is_some() {
            while let Some(&b) = dq.back() {
                let kb = keys[b].map(|x| x.0).unwrap_or(u64::MAX);
                let ki = keys[i].map(|x| x.0).unwrap_or(u64::MAX);
                // ties keep the leftmost (strict > pops)
                if kb > ki {
                    dq.pop_back();
                } else {
                    break;
                }
            }
            dq.push_back(i);
        }
        if i + 1 >= w {
            while let Some(&f) = dq.front() {
                if f + w <= i {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            if let Some(&f) = dq.front() {
                if let Some((_, kmer, strand)) = keys[f] {
                    // Adjacent windows selecting the same minimizer emit once.
                    if out.last().map(|m: &Minimizer| (m.kmer, m.pos)) != Some((kmer, f as u32)) {
                        out.push(Minimizer {
                            kmer,
                            pos: f as u32,
                            strand,
                        });
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_edge_kmers_are_bounded_canonical_deduplicated_and_query_only() {
        let seq = b"AGTCTGATCGTAGCTAGCATCGATCGTTAGCGATACGT";
        let params = SeedParams { k: 15, w: 5 };
        let mut scratch = SeedScratch::default();
        let mut out = Vec::new();
        minimizers_into(seq, params, &mut scratch, &mut out);
        let original = out.clone();
        let added = append_query_edge_kmers(seq.len(), params, &scratch, &mut out);
        assert_eq!(out.len(), original.len() + added);
        assert!(added <= 8);
        assert!(original.iter().all(|m| out.contains(m)));
        assert!(out.windows(2).all(|m| m[0].pos < m[1].pos));
        let n = seq.len() - 15 + 1;
        for pos in (0..4).chain(n - 4..n) {
            let mut expected = minimizers(&seq[pos..pos + 15], SeedParams { k: 15, w: 1 });
            expected[0].pos = pos as u32;
            assert!(out.contains(&expected[0]));
        }
        assert_eq!(
            append_query_edge_kmers(seq.len(), params, &scratch, &mut out),
            0
        );
        minimizers_into(seq, params, &mut scratch, &mut out);
        assert_eq!(
            out, original,
            "reference sketch and original extraction are unchanged"
        );
    }

    #[test]
    fn query_edge_kmers_reject_invalid_domains_and_skip_ambiguous_bases() {
        let seq = b"NNNTACGTTAGCTAGCTAGCATCGATCGATNGNNN";
        for params in [
            SeedParams { k: 15, w: 5 },
            SeedParams { k: 15, w: 1 },
            SeedParams { k: 15, w: 0 },
            SeedParams { k: 0, w: 5 },
            SeedParams { k: 33, w: 5 },
        ] {
            let mut scratch = SeedScratch::default();
            let mut out = Vec::new();
            // Populate first so early-return paths cannot reuse stale keys.
            minimizers_into(seq, SeedParams { k: 5, w: 2 }, &mut scratch, &mut out);
            minimizers_into(seq, params, &mut scratch, &mut out);
            let old = out.clone();
            let added = append_query_edge_kmers(seq.len(), params, &scratch, &mut out);
            if params.w <= 1 || params.k == 0 || params.k > 32 {
                assert_eq!(added, 0);
                assert_eq!(out, old);
            } else {
                for m in out {
                    assert!(seq[m.pos as usize..m.pos as usize + 15]
                        .iter()
                        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
                }
            }
        }
        for seq in [b"".as_slice(), b"ACGT", b"NNNNNNNNNNNNNNNNNNNN"] {
            let params = SeedParams { k: 15, w: 5 };
            let mut scratch = SeedScratch::default();
            let mut out = Vec::new();
            minimizers_into(seq, params, &mut scratch, &mut out);
            assert_eq!(
                append_query_edge_kmers(seq.len(), params, &scratch, &mut out),
                0
            );
            assert!(out.is_empty());
        }
    }

    #[test]
    fn query_edge_kmers_handle_overlapping_edges_and_reject_stale_lengths() {
        let seq = b"ACGTTACGTTACGTTACGTTACGTTACGTTACGTT";
        for k in [1, 15, 32] {
            let params = SeedParams { k, w: 50 };
            let mut scratch = SeedScratch::default();
            let mut out = Vec::new();
            minimizers_into(seq, params, &mut scratch, &mut out);
            append_query_edge_kmers(seq.len(), params, &scratch, &mut out);
            assert_eq!(out, minimizers(seq, SeedParams { k, w: 1 }));
            let before = out.clone();
            assert_eq!(
                append_query_edge_kmers(seq.len() + 1, params, &scratch, &mut out),
                0
            );
            assert_eq!(out, before);
        }
    }

    #[test]
    fn scratch_reuse_yields_independent_results() {
        let params = SeedParams { k: 5, w: 2 };
        let s1: &[u8] = b"ACGTACGTACGTTTACGGATC";
        let s2: &[u8] = b"TTTNGGGCCCCAAATGCA"; // N run: invalid windows
        let mut scratch = SeedScratch::default();
        let mut out = Vec::new();

        minimizers_into(s1, params, &mut scratch, &mut out);
        assert_eq!(out, minimizers(s1, params));
        assert!(!out.is_empty());

        // Second call through the same buffers: different input must give
        // exactly the fresh-call result (no state leak).
        minimizers_into(s2, params, &mut scratch, &mut out);
        assert_eq!(out, minimizers(s2, params));

        // Too-short input must clear the reused out buffer.
        minimizers_into(b"AC", params, &mut scratch, &mut out);
        assert!(out.is_empty());

        // And the first input again, reproducing the first result.
        minimizers_into(s1, params, &mut scratch, &mut out);
        assert_eq!(out, minimizers(s1, params));
    }
}
