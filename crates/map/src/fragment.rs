//! Observable paired read-through evidence, independent of reference or truth.

use crate::extend::CigarOp;
use crate::seed::Strand;
use serde::Serialize;

/// Restore the original query frame without adding scored/reference columns.
pub(crate) fn lift_cigar(cigar: &mut Vec<CigarOp>, strand: Strand, suffix: usize) {
    if suffix == 0 {
        return;
    }
    let n = u32::try_from(suffix).expect("raw read length fits CIGAR");
    let edge = if strand == Strand::Plus {
        cigar.last_mut()
    } else {
        cigar.first_mut()
    };
    if let Some(CigarOp::SoftClip(old)) = edge {
        *old += n;
    } else if strand == Strand::Plus {
        cigar.push(CigarOp::SoftClip(n));
    } else {
        cigar.insert(0, CigarOp::SoftClip(n));
    }
}

pub(crate) const ADAPTERS: [(&str, &[u8], &[u8]); 2] = [
    (
        "truseq",
        b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    ),
    ("nextera", b"CTGTCTCTTATACACATCT", b"CTGTCTCTTATACACATCT"),
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct FragmentEvidence {
    pub fragment_length: usize,
    pub overlap_bases: usize,
    pub viable_offsets: usize,
    pub adapter_family: &'static str,
    pub adapter_bases_r1: usize,
    pub adapter_bases_r2: usize,
}

/// Match QC overlap semantics, not alignment's ambiguity normalization.
/// Unknown bytes pass through unchanged; they may support competing offsets.
fn overlap_revcomp(seq: &[u8]) -> Vec<u8> {
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

pub(crate) fn infer(r1: &[u8], r2: &[u8]) -> Option<FragmentEvidence> {
    let max_t = r1.len().min(r2.len()).checked_sub(8)?;
    // Adapter support is a cheap pre-gate, not the overlap uniqueness test.
    let mut supported = Vec::new();
    for t in 31..=max_t {
        for &(family, a1, a2) in &ADAPTERS {
            let n1 = (r1.len() - t).min(a1.len());
            let n2 = (r2.len() - t).min(a2.len());
            if r1[t..t + n1] == a1[..n1] && r2[t..t + n2] == a2[..n2] {
                supported.push((t, family, n1, n2));
            }
        }
    }
    if supported.is_empty() {
        return None;
    }
    let rc = overlap_revcomp(r2);
    let mut viable = None;
    for offset in -(r2.len() as isize - 30) + 1..r1.len() as isize - 30 {
        let start1 = offset.max(0) as usize;
        let start2 = (-offset).max(0) as usize;
        let n = (r1.len() - start1).min(rc.len() - start2);
        let cap = 5.min(n / 10);
        let mut mismatches = 0;
        for i in 0..n {
            let a = r1[start1 + i];
            let b = rc[start2 + i];
            // Count every QC-viable alternative before imposing canonical
            // support on the sole winner. Matching unknown bytes, not just
            // Ns, can create competitors under QC's equality semantics.
            mismatches += usize::from(!a.eq_ignore_ascii_case(&b));
            if mismatches > cap {
                break;
            }
        }
        if mismatches <= cap {
            if viable.is_some() {
                return None; // Never choose one of several viable overlaps.
            }
            viable = Some((offset, n, mismatches));
        }
    }
    let (offset, n, mismatches) = viable?;
    if mismatches != 0 {
        return None;
    }
    let t = usize::try_from(offset + r2.len() as isize).ok()?;
    let &(_, family, n1, n2) = supported.iter().find(|x| x.0 == t)?;
    let mut counts = [0usize; 4];
    for base in &r1[..t] {
        let i = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        counts[i] += 1;
    }
    if counts.iter().filter(|&&c| c > 0).count() < 3 || counts.iter().any(|&c| c * 5 > t * 4) {
        return None;
    }
    Some(FragmentEvidence {
        fragment_length: t,
        overlap_bases: n,
        viable_offsets: 1,
        adapter_family: family,
        adapter_bases_r1: n1,
        adapter_bases_r2: n2,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn body(n: usize) -> Vec<u8> {
        let mut state = 0xe911da37u64;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    pub(crate) fn pair(fragment: &[u8], lengths: [usize; 2], family: usize) -> [Vec<u8>; 2] {
        let mut a = [fragment.to_vec(), ADAPTERS[family].1.to_vec()].concat();
        let mut b = [crate::align::revcomp(fragment), ADAPTERS[family].2.to_vec()].concat();
        a.resize(lengths[0], b'T');
        b.resize(lengths[1], b'T');
        [a, b]
    }

    #[test]
    fn paired_fragment_exact_supported_ends_and_unequal_lengths() {
        for n in [31, 50, 80, 110, 142] {
            for family in [0, 1] {
                for lengths in [[150, 150], [150, 170], [170, 150]] {
                    let [a, b] = pair(&body(n), lengths, family);
                    let evidence = infer(&a, &b).expect("supported read-through");
                    assert_eq!(evidence.fragment_length, n);
                    assert_eq!(evidence.overlap_bases, n);
                    assert_eq!(evidence.viable_offsets, 1);
                    assert_eq!(evidence.adapter_family, ADAPTERS[family].0);
                    assert!(evidence.adapter_bases_r1 >= 8 && evidence.adapter_bases_r2 >= 8);
                }
            }
        }
    }

    #[test]
    fn paired_fragment_rejects_missing_short_inconsistent_or_erroneous_evidence() {
        for n in [30, 143, 148, 150] {
            let [a, b] = pair(&body(n), [150, 150], 0);
            assert!(infer(&a, &b).is_none(), "unsupported length {n}");
        }
        let [a, b] = pair(&body(80), [150, 150], 0);
        for (mate, offset) in [(0, 79), (1, 79), (0, 84), (1, 84), (0, 80), (1, 80)] {
            let mut changed = [a.clone(), b.clone()];
            let base = &mut changed[mate][offset];
            *base = if *base == b'A' { b'C' } else { b'A' };
            assert!(infer(&changed[0], &changed[1]).is_none());
        }
        let mut unknown = a.clone();
        unknown[35] = b'N';
        assert!(infer(&unknown, &b).is_none());
        let mut absent = b.clone();
        absent[80..].fill(b'C');
        assert!(infer(&a, &absent).is_none());
        let mixed = pair(&body(80), [150, 150], 1);
        assert!(infer(&a, &mixed[1]).is_none());
        let other = pair(&body(90), [150, 150], 0);
        assert!(infer(&a, &other[1]).is_none());
    }

    #[test]
    fn paired_fragment_low_complexity_and_unsupported_layouts_abstain() {
        for body in [vec![b'A'; 80], b"AC".repeat(40), b"ACGT".repeat(20)] {
            let [a, b] = pair(&body, [150, 150], 0);
            assert!(infer(&a, &b).is_none());
        }
        let seq = body(220);
        assert!(infer(&seq[..150], &crate::align::revcomp(&seq[70..])).is_none());
        assert!(infer(&[], &[]).is_none());
        assert!(infer(&seq[..20], &seq[..20]).is_none());
    }

    #[test]
    fn paired_fragment_identical_observations_cannot_reveal_origin() {
        let [a, b] = pair(&body(80), [150, 150], 0);
        // These same bytes could instead be biological adapter-like sequence
        // in a different experimental interpretation. The inference cannot
        // distinguish the two origins; it has an explicit library assumption.
        let as_technical = infer(&a, &b);
        let as_biological = infer(&a.clone(), &b.clone());
        assert!(as_technical.is_some());
        assert_eq!(as_technical, as_biological);
    }

    #[test]
    fn paired_fragment_qc_ambiguous_competitor_must_not_be_hidden_by_n_filter() {
        let [mut a, mut b] = pair(&body(80), [150, 150], 0);
        a[80 + ADAPTERS[0].1.len()..].fill(b'N');
        b[80 + ADAPTERS[0].2.len()..].fill(b'N');
        // Besides the real -70 offset, QC accepts adapter-tail overlaps with
        // matching Ns. They must count against uniqueness even though they
        // cannot be the final canonical, exact supported fragment.
        assert!(infer(&a, &b).is_none());
    }

    #[test]
    fn paired_fragment_qc_unknown_bytes_must_not_hide_competing_offsets() {
        // QC preserves non-ACGT bytes, rather than interpreting their IUPAC
        // complement or normalizing them to N. Both cases must compete.
        for symbol in [b'R', b'r', b'Y', b'y', b'X', b'x'] {
            let [mut a, mut b] = pair(&body(80), [150, 150], 0);
            a[80 + ADAPTERS[0].1.len()..].fill(symbol);
            b[80 + ADAPTERS[0].2.len()..].fill(symbol);
            assert!(
                infer(&a, &b).is_none(),
                "QC-viable {} tails must count against uniqueness",
                char::from(symbol)
            );
        }
    }

    #[test]
    fn paired_fragment_counts_unsupported_positive_offsets_and_qc_error_cap() {
        for n in [30, 31, 35] {
            for errors in [0, 1, 3, 4] {
                let [mut a, mut b] = pair(&body(80), [150, 150], 0);
                let tail = body(n);
                a[150 - n..].copy_from_slice(&tail);
                b[150 - n..].copy_from_slice(&crate::align::revcomp(&tail));
                for i in 0..errors {
                    b[149 - i] = if b[149 - i] == b'A' { b'C' } else { b'A' };
                }
                // Unsupported offset 150-n still defeats uniqueness at
                // overlap 31/35 with up to three errors. Overlap 30 is outside
                // the QC search domain, and four errors exceed its cap.
                let competing = n > 30 && errors <= 3;
                assert_eq!(infer(&a, &b).is_none(), competing, "{n}/{errors}");
            }
        }
    }

    #[test]
    fn paired_fragment_raw_lift_merges_clips_and_preserves_reference_and_scores() {
        for strand in [Strand::Plus, Strand::Minus] {
            for cigar in [
                vec![CigarOp::Match(80)],
                vec![
                    CigarOp::SoftClip(2),
                    CigarOp::Match(75),
                    CigarOp::SoftClip(3),
                ],
                vec![
                    CigarOp::Match(40),
                    CigarOp::RefSkip(100),
                    CigarOp::Match(40),
                ],
                vec![CigarOp::Match(35), CigarOp::Ins(2), CigarOp::Match(43)],
            ] {
                let reference_span = |c: &[CigarOp]| {
                    c.iter()
                        .map(|op| match op {
                            CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => *n,
                            _ => 0,
                        })
                        .sum::<u32>()
                };
                let mut lifted = cigar.clone();
                lift_cigar(&mut lifted, strand, 70);
                assert_eq!(reference_span(&lifted), reference_span(&cigar));
                let query: u32 = lifted
                    .iter()
                    .map(|op| match op {
                        CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n,
                        _ => 0,
                    })
                    .sum();
                assert_eq!(query, 150);
                assert!(!lifted
                    .windows(2)
                    .any(|ops| matches!(ops, [CigarOp::SoftClip(_), CigarOp::SoftClip(_)])));
                assert!(
                    matches!(if strand == Strand::Plus { lifted.last() } else { lifted.first() },
                    Some(CigarOp::SoftClip(n)) if *n >= 70)
                );
            }
        }
    }

    #[test]
    fn paired_fragment_asymmetric_splice_clips_lift_only_the_sequenced_end() {
        let cigar = vec![
            CigarOp::SoftClip(3),
            CigarOp::Match(40),
            CigarOp::RefSkip(100),
            CigarOp::Match(35),
            CigarOp::SoftClip(2),
        ];
        let mut plus = cigar.clone();
        lift_cigar(&mut plus, Strand::Plus, 70);
        assert_eq!(
            plus,
            vec![
                CigarOp::SoftClip(3),
                CigarOp::Match(40),
                CigarOp::RefSkip(100),
                CigarOp::Match(35),
                CigarOp::SoftClip(72),
            ]
        );
        let mut minus = cigar;
        lift_cigar(&mut minus, Strand::Minus, 70);
        assert_eq!(
            minus,
            vec![
                CigarOp::SoftClip(73),
                CigarOp::Match(40),
                CigarOp::RefSkip(100),
                CigarOp::Match(35),
                CigarOp::SoftClip(2),
            ]
        );
    }
}
