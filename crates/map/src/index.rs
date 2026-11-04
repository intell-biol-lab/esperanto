//! Minimizer index over a reference genome.
//!
//! A read-only, shareable, zero-copy view. Legacy build/load APIs provide
//! process-lifetime storage; `OwnedIndex` views borrow scoped mmap storage.

use crate::fasta::{Base, Reference};
use crate::seed::{minimizers_from_codes, Minimizer, SeedParams, Strand};
use rayon::prelude::*;

/// Share fraction `HIGH_FREQ_FRAC` (parameter): a k-mer whose hit
/// count exceeds [`Index::freq_cutoff`] is treated as high-frequency.
pub const HIGH_FREQ_FRAC: f64 = 0.0002;

/// Position packing: `contig(31b) << 33 | strand(1b) << 32 | pos(32b)`.
#[inline]
pub fn pack_pos(contig: u32, strand: Strand, pos: u32) -> u64 {
    ((contig as u64) << 33) | ((strand as u64) << 32) | (pos as u64)
}

/// Inverse of [`pack_pos`]: `(contig, pos, strand)`.
#[inline]
pub fn unpack_pos(p: u64) -> (u32, u32, Strand) {
    let contig = (p >> 33) as u32;
    let pos = p as u32;
    let strand = if (p >> 32) & 1 == 1 {
        Strand::Minus
    } else {
        Strand::Plus
    };
    (contig, pos, strand)
}

impl Strand {
    /// XOR combination used to derive anchor strand from read/ref strands.
    #[inline]
    pub fn xor(a: Strand, b: Strand) -> Strand {
        if a == b {
            Strand::Plus
        } else {
            Strand::Minus
        }
    }
}

/// Read-only minimizer index.
pub struct Index<'a> {
    /// Seed parameters used at build time.
    pub params: SeedParams,
    /// Format version (paidx v1).
    pub version: u32,
    /// High-frequency cutoff: k-mers with `count > freq_cutoff` are skipped
    /// by anchor collection.
    pub freq_cutoff: u32,
    /// The reference these tables describe.
    pub reference: &'a Reference,
    /// Distinct canonical k-mers, ascending.
    pub kmers: &'a [u64],
    /// Start offset of each k-mer's position run in `positions`.
    pub offsets: &'a [u64],
    /// Length of each k-mer's position run.
    pub counts: &'a [u32],
    /// Packed positions, grouped by k-mer, ascending within each run.
    pub positions: &'a [u64],
}

/// Zero-allocation query result: a slice into the position table.
pub struct IndexHit<'a> {
    /// Number of reference occurrences of the queried k-mer.
    pub count: u32,
    /// Packed positions (ascending).
    pub positions: &'a [u64],
}

impl Index<'static> {
    /// Build an index from a parsed reference (per-contig rayon parallel).
    ///
    /// This legacy constructor retains process-lifetime tables/reference.
    pub fn build(reference: Reference, params: SeedParams) -> Index<'static> {
        let k = params.k;
        let w = params.w;

        let per_contig: Vec<Vec<(u64, u64)>> = reference
            .contigs
            .par_iter()
            .enumerate()
            .map(|(ci, contig)| {
                let codes = contig_codes(contig);
                let mins = minimizers_from_codes(&codes, k, w);
                mins.into_iter()
                    .map(|m| (m.kmer, pack_pos(ci as u32, m.strand, m.pos)))
                    .collect()
            })
            .collect();

        let mut all: Vec<(u64, u64)> = per_contig.into_iter().flatten().collect();
        all.sort_unstable();
        all.dedup();

        let mut kmers: Vec<u64> = Vec::new();
        let mut offsets: Vec<u64> = Vec::new();
        let mut counts: Vec<u32> = Vec::new();
        let mut positions: Vec<u64> = Vec::new();
        let mut i = 0usize;
        while i < all.len() {
            let kmer = all[i].0;
            let start = i;
            while i < all.len() && all[i].0 == kmer {
                positions.push(all[i].1);
                i += 1;
            }
            kmers.push(kmer);
            offsets.push(start as u64);
            counts.push((i - start) as u32);
        }

        let freq_cutoff = compute_cutoff(&counts);

        let reference: &'static Reference = Box::leak(Box::new(reference));
        let kmers: &'static [u64] = Box::leak(kmers.into_boxed_slice());
        let offsets: &'static [u64] = Box::leak(offsets.into_boxed_slice());
        let counts: &'static [u32] = Box::leak(counts.into_boxed_slice());
        let positions: &'static [u64] = Box::leak(positions.into_boxed_slice());

        Index {
            params,
            version: 1,
            freq_cutoff,
            reference,
            kmers,
            offsets,
            counts,
            positions,
        }
    }
}

impl<'a> Index<'a> {
    /// Query a canonical k-mer; binary search, zero allocation.
    pub fn query(&self, kmer: u64) -> IndexHit<'a> {
        let pos = self.kmers.binary_search(&kmer);
        match pos {
            Ok(i) => {
                let start = self.offsets[i] as usize;
                let count = self.counts[i] as usize;
                IndexHit {
                    count: self.counts[i],
                    positions: &self.positions[start..start + count],
                }
            }
            Err(_) => IndexHit {
                count: 0,
                positions: &[],
            },
        }
    }

    /// Whether `count` marks a k-mer as high-frequency (skip in anchoring).
    pub fn is_high_freq(&self, count: u32) -> bool {
        count > self.freq_cutoff
    }
}

/// Frequency cutoff from per-k-mer counts: sort descending, take the value at
/// rank `len * HIGH_FREQ_FRAC` (capped at `len - 1`), minimum 1. Empty input
/// yields `u32::MAX` (nothing is high-frequency).
pub fn compute_cutoff(counts: &[u32]) -> u32 {
    if counts.is_empty() {
        return u32::MAX;
    }
    let mut sorted: Vec<u32> = counts.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let rank = ((sorted.len() as f64) * HIGH_FREQ_FRAC) as usize;
    let rank = rank.min(sorted.len() - 1);
    sorted[rank].max(1)
}

fn contig_codes(contig: &crate::fasta::Contig) -> Vec<u64> {
    const N_CODE: u64 = 0xFF;
    contig
        .slice_ascii(0, contig.len)
        .into_iter()
        .map(|b| match Base::from_ascii(b).code() {
            Some(c) => c as u64,
            None => N_CODE,
        })
        .collect()
}

/// qpos for an anchor, keyed on the ANCHOR strand (not the minimizer's):
/// Plus → `m.pos`, Minus → `read_len − m.pos − k`.
#[inline]
fn anchor_qpos(astrand: Strand, m: &Minimizer, read_len: u32, k: u32) -> u32 {
    if astrand == Strand::Plus {
        m.pos
    } else {
        read_len - m.pos - k
    }
}

/// Collect anchors for plain read minimizers.
///
/// Skips k-mers with `count == 0` or `count > occ_cap` (all positions of a
/// passing k-mer are emitted). Anchor strand is the XOR of the read and
/// reference minimizer strands (same → Plus, different → Minus), and qpos is
/// derived from the anchor strand: Plus anchor → `m.pos`, Minus anchor →
/// `read_len − m.pos − k`. Returns the anchors and the number of read
/// minimizers that passed the cap.
pub fn collect_anchors(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
) -> (Vec<crate::seed::Anchor>, usize) {
    let mut anchors = Vec::new();
    let seeds_hit = collect_anchors_into(index, read_mins, read_len, k, occ_cap, &mut anchors);
    (anchors, seeds_hit)
}

/** [`collect_anchors`] into a reusable buffer (`anchors` is cleared, then
 * filled; same anchors, same order).
 */
pub fn collect_anchors_into(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
    anchors: &mut Vec<crate::seed::Anchor>,
) -> usize {
    anchors.clear();
    let mut seeds_hit = 0usize;
    for m in read_mins {
        let hit = index.query(m.kmer);
        if hit.count == 0 || hit.count > occ_cap {
            continue;
        }
        seeds_hit += 1;
        for &p in hit.positions {
            let (contig, rpos, rstrand) = unpack_pos(p);
            let astrand = Strand::xor(m.strand, rstrand);
            anchors.push(crate::seed::Anchor {
                qpos: anchor_qpos(astrand, m, read_len, k),
                rpos,
                contig,
                strand: astrand,
            });
        }
    }
    seeds_hit
}

/// Same-domain background evidence for additional ordinary context proposals.
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct ContextBackgroundBudget {
    pub score: i32,
    pub rejected_candidates: usize,
    pub terminal_prefix_bases: usize,
    pub terminal_suffix_bases: usize,
    pub rejected_terminal_candidates: usize,
}

/// Additional context verification, not ordinary retained anchors.
/// Its presence denotes a conservative, ungapped-only limited search.
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct RepeatContextBudget {
    pub query_length: usize,
    pub eligible_seed_lists: usize,
    pub scanned_seed_lists: usize,
    pub omitted_seed_lists: usize,
    pub eligible_postings: usize,
    pub scanned_postings: usize,
    pub unique_intervals: usize,
    pub positive_candidates: usize,
    pub retained_candidates: usize,
    pub omitted_candidates: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub background: Option<ContextBackgroundBudget>,
}

#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct AnchorBudget {
    pub eligible_anchors: usize,
    pub retained_anchors: usize,
    pub omitted_seeds: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub repeat_context: Option<RepeatContextBudget>,
}

impl AnchorBudget {
    pub fn limited(self) -> bool {
        self.omitted_seeds != 0 || self.repeat_context.is_some()
    }
}

/// Collect whole seed hit lists within `budget`, preferring rarer seeds only
/// on overflow. Every locus for an admitted seed is retained; no coordinate
/// prefix sampling is allowed. The caller must not infer uniqueness from a
/// limited search. Complete inputs preserve the original collector's order.
#[allow(clippy::too_many_arguments)]
pub fn collect_anchors_budgeted_into(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
    budget: usize,
    anchors: &mut Vec<crate::seed::Anchor>,
) -> AnchorBudget {
    anchors.clear();
    let mut seeds = Vec::new();
    let mut eligible_anchors = 0usize;
    for (i, m) in read_mins.iter().enumerate() {
        let count = index.query(m.kmer).count;
        if count > 0 && count <= occ_cap {
            seeds.push((i, count as usize));
            eligible_anchors = eligible_anchors.saturating_add(count as usize);
        }
    }
    if eligible_anchors <= budget {
        collect_anchors_into(index, read_mins, read_len, k, occ_cap, anchors);
        return AnchorBudget {
            eligible_anchors,
            retained_anchors: anchors.len(),
            omitted_seeds: 0,
            repeat_context: None,
        };
    }
    seeds.sort_unstable_by_key(|&(i, count)| {
        let m = read_mins[i];
        (count, m.pos, m.kmer, m.strand as u8)
    });
    let mut omitted_seeds = 0;
    for (i, count) in seeds {
        if count > budget - anchors.len() {
            omitted_seeds += 1;
            continue;
        }
        let m = &read_mins[i];
        for &p in index.query(m.kmer).positions {
            let (contig, rpos, rstrand) = unpack_pos(p);
            let strand = Strand::xor(m.strand, rstrand);
            anchors.push(crate::seed::Anchor {
                qpos: anchor_qpos(strand, m, read_len, k),
                rpos,
                contig,
                strand,
            });
        }
    }
    AnchorBudget {
        eligible_anchors,
        retained_anchors: anchors.len(),
        omitted_seeds,
        repeat_context: None,
    }
}

/// Reverse-complement a 2-bit k-mer of length `k`.
fn revcomp_kmer(kmer: u64, k: u32) -> u64 {
    let mut out: u64 = 0;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 3;
        out |= (3 - base) << (2 * (k - 1 - i));
    }
    out
}

/// Canonical form (min of k-mer and revcomp) with orientation flag.
fn canonical(kmer: u64, k: u32) -> (u64, Strand) {
    let rc = revcomp_kmer(kmer, k);
    if kmer <= rc {
        (kmer, Strand::Plus)
    } else {
        (rc, Strand::Minus)
    }
}

/// Emit anchors for one hitting variant round.
///
/// `take` caps the positions consumed;
/// `rev` selects the round: fwd round anchors carry the reference strand,
/// rev round anchors carry its flip. qpos always follows the anchor strand.
fn emit_variant_round(
    anchors: &mut Vec<crate::seed::Anchor>,
    positions: &[u64],
    m: &Minimizer,
    read_len: u32,
    k: u32,
    take: usize,
    rev: bool,
) {
    for &p in positions.iter().take(take) {
        let (contig, rpos, rstrand) = unpack_pos(p);
        let astrand = if rev { rstrand.flip() } else { rstrand };
        anchors.push(crate::seed::Anchor {
            qpos: anchor_qpos(astrand, m, read_len, k),
            rpos,
            contig,
            strand: astrand,
        });
    }
}

/// Single-site G→A / C→T edit variants (editing-aware anchor collection).
///
/// For each read minimizer, the k-mer is packed in fwd orientation (site i
/// lives at 2-bit slot `k − 1 − i`); every editable site runs TWO rounds —
/// fwd-orientation variant (want = alt) and rev-orientation variant (mirrored
/// site, want = 3 − alt). Both rounds canonicalize to the same key and BOTH
/// query and BOTH emit: a hitting site contributes its first 64 positions as
/// anchors with the reference strand (fwd round), then the same first 64 as
/// anchors with the flipped strand (rev round), so every hit position yields
/// exactly one `(m.pos, Plus)` and one `(read_len − m.pos − k, Minus)` anchor
/// Single-site G→A / C→T edit variants (editing-aware anchor collection).
///
/// Per editable site i (read-orientation k-mer base G/C → alt A/T), TWO
/// rounds run, each producing a DIFFERENT variant key:
/// - fwd round: edit slot (k−1−i) of the read-forward packing ← alt (edits
///   read position i); anchors carry the reference strand;
/// - rev round: edit slot (k−1−i) of the rev packing ← 3−alt (edits the
///   MIRROR read position k−1−i); skipped when the rev slot already holds
///   3−alt (no-op: mirror read base already equals alt); anchors carry the
///   flipped reference strand.
///
/// qpos follows the anchor strand (Plus → m.pos, Minus → read_len−m.pos−k).
/// Each hitting round takes the first 64 positions and increments var_hits.
pub fn collect_anchors_edit_variants(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
) -> (Vec<crate::seed::Anchor>, usize) {
    let mut anchors = Vec::new();
    let mut var_hits = 0usize;
    for m in read_mins {
        // k-mer as it reads on the read's forward strand
        let fwd = if m.strand == Strand::Minus {
            revcomp_kmer(m.kmer, k)
        } else {
            m.kmer
        };
        let rev = revcomp_kmer(fwd, k);
        for i in 0..k as usize {
            let slot = (k as usize - 1) - i;
            let base = (fwd >> (2 * slot)) & 3;
            let alt = match base {
                2 => 0u64, // G -> A
                1 => 3u64, // C -> T
                _ => continue,
            };
            // fwd round
            let v = (fwd & !(3 << (2 * slot))) | (alt << (2 * slot));
            let key = canonical(v, k).0;
            let hit = index.query(key);
            if hit.count > 0 && hit.count <= occ_cap {
                var_hits += 1;
                emit_variant_round(&mut anchors, hit.positions, m, read_len, k, 64, false);
            }
            // rev round (mirror-position edit; skip no-ops)
            let want = 3 - alt;
            let cur = (rev >> (2 * slot)) & 3;
            if cur == want {
                continue;
            }
            let v = (rev & !(3 << (2 * slot))) | (want << (2 * slot));
            let key = canonical(v, k).0;
            let hit = index.query(key);
            if hit.count > 0 && hit.count <= occ_cap {
                var_hits += 1;
                emit_variant_round(&mut anchors, hit.positions, m, read_len, k, 64, true);
            }
        }
    }
    (anchors, var_hits)
}

#[cfg(test)]
mod budget_tests {
    use super::*;
    use crate::seed::Anchor;

    fn fixture() -> (Index<'static>, Vec<Minimizer>) {
        let reference =
            crate::fasta::parse_fasta_bytes(format!(">chrT\n{}\n", "ACGT".repeat(500)).as_bytes())
                .unwrap();
        let counts = vec![4, 2, 3, 6];
        let mut positions = Vec::new();
        let mut offsets = Vec::new();
        for (i, count) in counts.iter().enumerate() {
            offsets.push(positions.len() as u64);
            for j in 0..*count {
                positions.push(pack_pos(
                    0,
                    if j % 2 == 0 {
                        Strand::Plus
                    } else {
                        Strand::Minus
                    },
                    100 * i as u32 + 10 * j,
                ));
            }
        }
        let index = Index {
            params: SeedParams { k: 15, w: 5 },
            version: 1,
            freq_cutoff: 1,
            reference: Box::leak(Box::new(reference)),
            kmers: Box::leak(vec![1, 2, 3, 4].into_boxed_slice()),
            counts: Box::leak(counts.into_boxed_slice()),
            offsets: Box::leak(offsets.into_boxed_slice()),
            positions: Box::leak(positions.into_boxed_slice()),
        };
        let mins = [
            (1, 10, Strand::Plus),
            (3, 50, Strand::Minus),
            (2, 30, Strand::Plus),
            (4, 70, Strand::Plus),
        ]
        .map(|(kmer, pos, strand)| Minimizer { kmer, pos, strand })
        .to_vec();
        (index, mins)
    }

    #[test]
    fn budget_complete_input_preserves_anchor_order_and_strands() {
        let (index, mins) = fixture();
        let (expected, _) = collect_anchors(&index, &mins, 150, 15, 5);
        for budget in [9, 10, 20000] {
            let mut actual = Vec::new();
            let status =
                collect_anchors_budgeted_into(&index, &mins, 150, 15, 5, budget, &mut actual);
            assert_eq!(actual, expected);
            assert!(!status.limited());
            assert_eq!(status.eligible_anchors, 9);
            assert_eq!(status.omitted_seeds, 0);
        }
    }

    #[test]
    fn budget_prefers_rare_seeds_and_never_samples_reference_coordinates() {
        let (index, mins) = fixture();
        for budget in [5, 8] {
            let mut actual = Vec::<Anchor>::new();
            let status =
                collect_anchors_budgeted_into(&index, &mins, 150, 15, 5, budget, &mut actual);
            let expected_mins = [mins[2], mins[1]];
            assert_eq!(
                actual,
                collect_anchors(&index, &expected_mins, 150, 15, 5).0
            );
            assert_eq!(actual.len(), 5);
            assert!(actual.len() <= budget);
            assert!(status.limited());
            assert_eq!(status.omitted_seeds, 1);
            assert_eq!(status.eligible_anchors, 9);
        }
    }

    #[test]
    fn budget_zero_and_too_small_abstain_without_partial_seed_lists() {
        let (index, mins) = fixture();
        for budget in [0, 1] {
            let mut actual = vec![Anchor {
                contig: 0,
                qpos: 0,
                rpos: 0,
                strand: Strand::Plus,
            }];
            let status =
                collect_anchors_budgeted_into(&index, &mins, 150, 15, 5, budget, &mut actual);
            assert!(actual.is_empty());
            assert_eq!(status.omitted_seeds, 3);
            assert!(status.limited());
        }
    }

    #[test]
    fn budget_selection_is_deterministic_when_seed_input_is_permuted() {
        let (index, mut mins) = fixture();
        let mut expected = Vec::new();
        collect_anchors_budgeted_into(&index, &mins, 150, 15, 5, 5, &mut expected);
        mins.reverse();
        let mut actual = Vec::new();
        collect_anchors_budgeted_into(&index, &mins, 150, 15, 5, 5, &mut actual);
        assert_eq!(actual, expected);
    }
}
