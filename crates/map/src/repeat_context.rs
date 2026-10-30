//! Bounded ordinary ungapped hypotheses from complete high-frequency lists.
//! This adds candidates, never replaces old splice/indel candidates or DP.
use std::collections::BTreeSet;

use crate::align::{revcomp, MID_OCC_CAP, REFUSE_MIN_ALIGNED_FRAC};
use crate::extend::{exempt_dir, CigarOp, ExtendParams};
use crate::index::{unpack_pos, ContextBackgroundBudget, Index, RepeatContextBudget};
use crate::seed::{minimizers, SeedParams, Strand};

pub(crate) const POSTING_BUDGET: usize = 20_000;
pub(crate) const PROPOSAL_BUDGET: usize = 8;
pub(crate) const MAX_QUERY: usize = 256;

pub(crate) fn generate_budgeted(
    index: &Index,
    query: &[u8],
    ep: ExtendParams,
    background: bool,
    postings: usize,
    proposals: usize,
) -> Option<(Vec<Proposal>, RepeatContextBudget)> {
    generate_impl(index, query, ep, background, postings, proposals)
}

pub(crate) struct Proposal {
    pub contig: u32,
    pub pos: u32,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub score: i32,
    ordinary: u32,
    informative: u32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Local {
    start: usize,
    end: usize,
    score: i32,
    ordinary: u32,
    informative: u32,
}

fn local(query: &[u8], reference: &[u8], ep: ExtendParams) -> Local {
    debug_assert_eq!(query.len(), reference.len());
    let mut best = Local::default();
    let mut cur = Local::default();
    for (i, (&q, &r)) in query.iter().zip(reference).enumerate() {
        // Preserve a zero-score informative editing block, but never buy
        // coverage by crossing a zero-score ordinary or ambiguous prefix.
        if cur.score < 0 || (cur.score == 0 && (cur.ordinary > 0 || cur.informative == 0)) {
            cur = Local {
                start: i,
                ..Local::default()
            };
        }
        cur.end = i + 1;
        if q != b'N' && r != b'N' {
            cur.informative += 1;
            if q == r {
                cur.score += ep.match_score;
            } else {
                cur.score += ep.mismatch;
                cur.ordinary += u32::from(exempt_dir(r, q).is_none());
            }
        }
        if cur.score > best.score
            || (cur.score == best.score
                && (cur.ordinary < best.ordinary
                    || (cur.ordinary == best.ordinary && cur.informative > best.informative)))
        {
            best = cur;
        }
    }
    best
}

/// Best strict local score against a coordinate-free constant-base reference.
/// This is a competing explanation, not a technical-origin classifier.
fn background_score(query: &[u8], ep: ExtendParams) -> i32 {
    b"ACGT"
        .iter()
        .map(|&base| {
            let mut best = 0;
            let mut score = 0;
            for &q in query {
                let s = if q == b'N' {
                    0
                } else if q == base {
                    ep.match_score
                } else {
                    ep.mismatch
                };
                score = (score + s).max(0);
                best = best.max(score);
            }
            best
        })
        .max()
        .unwrap_or(0)
}

/// Positive local M/S proposals from a bounded union of complete hit lists.
/// `None` means no additional search; `Some` is conservatively limited even
/// when every selected list fits, because this family verifies ungapped
/// diagonals only and does not establish a calibrated complete search.
#[cfg(test)]
pub(crate) fn generate(
    index: &Index,
    query: &[u8],
    ep: ExtendParams,
) -> Option<(Vec<Proposal>, RepeatContextBudget)> {
    generate_with_background(index, query, ep, true)
}

#[cfg(test)]
pub(crate) fn generate_with_background(
    index: &Index,
    query: &[u8],
    ep: ExtendParams,
    background: bool,
) -> Option<(Vec<Proposal>, RepeatContextBudget)> {
    generate_impl(
        index,
        query,
        ep,
        background,
        POSTING_BUDGET,
        PROPOSAL_BUDGET,
    )
}

fn generate_impl(
    index: &Index,
    query: &[u8],
    ep: ExtendParams,
    background: bool,
    posting_budget: usize,
    proposal_budget: usize,
) -> Option<(Vec<Proposal>, RepeatContextBudget)> {
    let k = index.params.k as usize;
    if k == 0 || k > 32 || query.len() < k || query.len() > MAX_QUERY {
        return None;
    }
    let mut seeds: Vec<_> = minimizers(query, SeedParams { k: k as u32, w: 1 })
        .into_iter()
        .filter_map(|m| {
            let count = index.query(m.kmer).count as usize;
            (count > MID_OCC_CAP as usize).then_some((count, m))
        })
        .collect();
    if seeds.is_empty() {
        return None;
    }
    seeds.sort_unstable_by_key(|(count, m)| (*count, m.pos, m.kmer, m.strand as u8));
    let mut budget = RepeatContextBudget {
        query_length: query.len(),
        eligible_seed_lists: seeds.len(),
        eligible_postings: seeds.iter().map(|(n, _)| n).sum(),
        background: background.then(|| ContextBackgroundBudget {
            score: background_score(query, ep),
            rejected_candidates: 0,
        }),
        ..RepeatContextBudget::default()
    };
    let rc = revcomp(query);
    let mut seen = BTreeSet::new();
    let mut emitted = BTreeSet::new();
    let mut proposals = Vec::new();
    for (count, m) in seeds {
        if count > posting_budget - budget.scanned_postings {
            budget.omitted_seed_lists += 1;
            continue;
        }
        budget.scanned_seed_lists += 1;
        budget.scanned_postings += count;
        for &packed in index.query(m.kmer).positions {
            let (ci, rp, rs) = unpack_pos(packed);
            let strand = Strand::xor(m.strand, rs);
            let qp = if strand == Strand::Plus {
                m.pos as usize
            } else {
                query.len() - m.pos as usize - k
            };
            let Some(start) = rp.checked_sub(qp as u32) else {
                continue;
            };
            let contig = &index.reference.contigs[ci as usize];
            let Some(end) = start.checked_add(query.len() as u32) else {
                continue;
            };
            let key = (ci, start, strand as u8);
            if end > contig.len || emitted.contains(&key) {
                continue;
            }
            // A first seed outside the best local interval does not qualify
            // the proposal, but cannot erase a later inside seed on the same
            // diagonal. Unique interval accounting remains truthful.
            budget.unique_intervals += usize::from(seen.insert(key));
            let q = if strand == Strand::Plus { query } else { &rc };
            let reference = contig.slice_ascii(start, end);
            // Do not turn normalization/canonicalization into fictitious
            // byte-exact anchor support.
            if q[qp..qp + k] != reference[qp..qp + k] {
                continue;
            }
            let best = local(q, &reference, ep);
            if best.score <= 0 || best.start > qp || best.end < qp + k {
                continue;
            }
            emitted.insert(key);
            if let Some(background) = &mut budget.background {
                if best.score <= background.score {
                    background.rejected_candidates += 1;
                    continue;
                }
            }
            budget.positive_candidates += 1;
            let mut cigar = Vec::with_capacity(3);
            if best.start > 0 {
                cigar.push(CigarOp::SoftClip(best.start as u32));
            }
            cigar.push(CigarOp::Match((best.end - best.start) as u32));
            if best.end < query.len() {
                cigar.push(CigarOp::SoftClip((query.len() - best.end) as u32));
            }
            proposals.push(Proposal {
                contig: ci,
                pos: start + best.start as u32,
                strand,
                cigar,
                score: best.score,
                ordinary: best.ordinary,
                informative: best.informative,
            });
        }
    }
    // Preserve a primary-eligible fallback even when short local competitors
    // have higher scores. All omitted competitors remain limited evidence.
    let eligible =
        |p: &Proposal| f64::from(p.informative) >= REFUSE_MIN_ALIGNED_FRAC * query.len() as f64;
    proposals.sort_unstable_by(|a, b| {
        eligible(b)
            .cmp(&eligible(a))
            .then(b.score.cmp(&a.score))
            .then(a.ordinary.cmp(&b.ordinary))
            .then(b.informative.cmp(&a.informative))
            .then(a.contig.cmp(&b.contig))
            .then(a.pos.cmp(&b.pos))
            .then((a.strand as u8).cmp(&(b.strand as u8)))
    });
    proposals.truncate(proposal_budget);
    budget.retained_candidates = proposals.len();
    budget.omitted_candidates = budget.positive_candidates - proposals.len();
    Some((proposals, budget))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_scan_and_retention_budgets_preserve_whole_lists() {
        let query = b"ACGTTGCAAGTCGAT";
        let reference = [vec![b'N'; 32], query.to_vec(), vec![b'N'; 32]]
            .concat()
            .repeat(6000);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">repeat\n".as_slice(), &reference, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 1 },
        );
        let ep = ExtendParams::default();
        let (none, skipped) = generate_budgeted(&index, query, ep, true, 5999, 8).unwrap();
        assert!(
            none.is_empty(),
            "never sample a coordinate prefix of a list"
        );
        assert_eq!(skipped.scanned_postings, 0);
        assert_eq!(skipped.omitted_seed_lists, 1);
        let (small, a) = generate_budgeted(&index, query, ep, true, 6000, 8).unwrap();
        let (large, b) = generate_budgeted(&index, query, ep, true, 6000, 32).unwrap();
        assert_eq!((small.len(), large.len()), (8, 32));
        assert_eq!((a.scanned_postings, b.scanned_postings), (6000, 6000));
        assert_eq!(a.positive_candidates, b.positive_candidates);
        assert_eq!(a.omitted_candidates - b.omitted_candidates, 24);
        let (defaults, _) = generate(&index, query, ep).unwrap();
        assert_eq!(
            defaults
                .iter()
                .map(|p| (p.contig, p.pos, p.score))
                .collect::<Vec<_>>(),
            small
                .iter()
                .map(|p| (p.contig, p.pos, p.score))
                .collect::<Vec<_>>()
        );
        // Exhausting these ungapped lists does not make the evidence complete.
        assert!(crate::index::AnchorBudget {
            repeat_context: Some(b),
            ..Default::default()
        }
        .limited());
    }

    #[test]
    fn context_background_score_matches_exhaustive_subinterval_oracle() {
        let mut state = 0x175831_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        assert_eq!(background_score(b"", ExtendParams::default()), 0);
        for i in 0..2500 {
            let n = (next() % 48 + 1) as usize;
            let query: Vec<_> = (0..n).map(|_| b"ACGTN"[(next() % 5) as usize]).collect();
            let ep = if i % 2 == 0 {
                ExtendParams::default()
            } else {
                ExtendParams {
                    match_score: 3,
                    mismatch: -2,
                    ..ExtendParams::default()
                }
            };
            let mut expected = 0;
            for &base in b"ACGT" {
                for start in 0..n {
                    let mut score = 0;
                    for &q in &query[start..] {
                        score += if q == b'N' {
                            0
                        } else if q == base {
                            ep.match_score
                        } else {
                            ep.mismatch
                        };
                        expected = expected.max(score);
                    }
                }
            }
            assert_eq!(background_score(&query, ep), expected);
            assert_eq!(background_score(&revcomp(&query), ep), expected);
        }
        // The best constant explanation can bridge a mismatch; it is not just
        // a longest-run heuristic. N contributes neither score nor information.
        assert_eq!(background_score(b"AAAACAAAA", ExtendParams::default()), 12);
        assert_eq!(background_score(b"AAANAAA", ExtendParams::default()), 12);
    }

    #[test]
    fn context_background_equality_rejects_and_counts_once_per_diagonal() {
        let query = vec![b'A'; 15];
        // Native sketch: one valid kmer per separated copy, no forced cap.
        let repeated = [vec![b'N'; 32], query.clone(), vec![b'N'; 32]]
            .concat()
            .repeat(6000);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">repeat\n".as_slice(), &repeated, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let (off, b0) =
            generate_with_background(&index, &query, ExtendParams::default(), false).unwrap();
        assert!(!off.is_empty());
        assert_eq!(off[0].score, 30);
        assert!(b0.background.is_none());
        let (on, b1) = generate(&index, &query, ExtendParams::default()).unwrap();
        assert!(on.is_empty(), "a constant-base tie provides no specificity");
        let bg = b1.background.unwrap();
        assert_eq!(bg.score, 30);
        assert_eq!(bg.rejected_candidates, b0.positive_candidates);
        assert_eq!(b1.positive_candidates, 0);
        assert_eq!(b1.scanned_postings, b0.scanned_postings);
        assert_eq!(b1.unique_intervals, b0.unique_intervals);
    }

    #[test]
    fn local_context_ties_keep_neutral_edits_but_not_ordinary_or_ambiguous_tails() {
        let ep = ExtendParams::default();
        let r = b"CCAACGTACC";
        let edited = local(b"CCGACGTGCC", r, ep);
        assert_eq!(
            (edited.start, edited.end, edited.score, edited.ordinary),
            (0, 10, 8, 0)
        );
        let ordinary = local(b"CCTACGTTCC", r, ep);
        assert_eq!(
            (
                ordinary.start,
                ordinary.end,
                ordinary.score,
                ordinary.ordinary
            ),
            (3, 7, 8, 0)
        );
        let ambiguous = local(b"NNACGTNN", b"NNACGTNN", ep);
        assert_eq!(
            (ambiguous.start, ambiguous.end, ambiguous.informative),
            (2, 6, 4)
        );
        let negative = local(b"ACGTCCC", b"ACGTAAA", ep);
        assert_eq!((negative.start, negative.end, negative.score), (0, 4, 8));
    }

    #[test]
    fn positive_technical_match_is_an_explicit_limit_not_a_promised_discriminator() {
        // Identical score evidence cannot reveal biological vs technical origin.
        let p = local(b"ACGTCAAA", b"ACGTAAAA", ExtendParams::default());
        assert_eq!((p.start, p.end, p.score, p.ordinary), (0, 8, 10, 1));
    }

    #[test]
    fn query_bounds_and_no_high_frequency_domain_do_not_create_search_evidence() {
        let (index, seq) = crate::align::repeat_context_tests::fixture(false, false, 6000);
        assert!(generate(&index, &seq[..14], ExtendParams::default()).is_none());
        assert!(generate(&index, &vec![b'N'; 257], ExtendParams::default()).is_none());
        assert!(generate(&index, &vec![b'N'; 112], ExtendParams::default()).is_none());
        let bg = index.reference.contigs[1].slice_ascii(200, 350);
        assert!(generate(&index, &bg, ExtendParams::default()).is_none());
    }

    #[test]
    fn local_context_matches_bruteforce_common_order_on_generated_diagonals() {
        let mut state = 0x7517_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..3000 {
            let n = (next() % 24 + 1) as usize;
            let mut q = Vec::new();
            let mut r = Vec::new();
            for _ in 0..n {
                let a = b"ACGTN"[(next() % 5) as usize];
                r.push(a);
                q.push(if next() & 1 == 0 {
                    a
                } else {
                    b"ACGTN"[(next() % 5) as usize]
                });
            }
            let actual = local(&q, &r, ExtendParams::default());
            let mut best = (0, std::cmp::Reverse(0_u32), 0_u32);
            for s in 0..n {
                let (mut score, mut mm, mut covered) = (0, 0_u32, 0_u32);
                for e in s..n {
                    if q[e] != b'N' && r[e] != b'N' {
                        covered += 1;
                        if q[e] == r[e] {
                            score += 2;
                        } else {
                            score -= 4;
                            mm += u32::from(exempt_dir(r[e], q[e]).is_none());
                        }
                    }
                    best = best.max((score, std::cmp::Reverse(mm), covered));
                }
            }
            assert_eq!(actual.score, best.0);
            if actual.score > 0 {
                assert_eq!(
                    (
                        actual.score,
                        std::cmp::Reverse(actual.ordinary),
                        actual.informative
                    ),
                    best,
                    "query={q:?} reference={r:?}"
                );
            }
        }
    }

    #[test]
    fn later_seed_inside_local_interval_is_not_erased_by_an_earlier_outside_seed() {
        let mut state = 149_u64;
        let body: Vec<_> = (0..112)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let mut query = body.clone();
        for p in 15..31 {
            query[p] = b"ACGT"[(b"ACGT".iter().position(|&b| b == query[p]).unwrap() + 1) % 4];
        }
        let reference = [vec![b'N'; 32], body, vec![b'N'; 32]].concat().repeat(6000);
        let index = Index::build(
            crate::fasta::parse_fasta_bytes(&[b">repeat\n".as_slice(), &reference, b"\n"].concat())
                .unwrap(),
            SeedParams { k: 15, w: 5 },
        );
        let (p, b) = generate(&index, &query, ExtendParams::default()).unwrap();
        assert_eq!(b.scanned_seed_lists, 3);
        assert!(
            !p.is_empty(),
            "an outside first seed cannot hide an inside later seed"
        );
        assert!(p
            .iter()
            .all(|p| p.cigar.first() == Some(&CigarOp::SoftClip(31))));
    }
}
