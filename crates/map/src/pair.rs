//! Paired-end relationship and insert-size statistics.

use crate::extend::CigarOp;
use crate::mapq::ReadAlignment;

/// Maximum reference span for a proper pair ( 1000).
pub const PROPER_SPAN_MAX: u32 = 1000;

/// Relationship of two mates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairRelation {
    /// Mates on opposite strands (FR or RF) on the same contig within the
    /// transcript-aware span cap.
    Proper,
    /// Anything else.
    Discordant,
}

/// Classify a mate pair: proper iff same contig, mates on opposite strands
/// (FR or RF — library strand specificity flips which orientation dominates,
/// but either is a valid inward-facing pair), and the transcript-aware span
/// (genomic span minus intronic N-skips) is within the span cap.
pub fn relate(r1: &ReadAlignment, r2: &ReadAlignment) -> PairRelation {
    if r1.contig != r2.contig {
        return PairRelation::Discordant;
    }
    // Same-strand (FF/RR) is discordant; opposite-strand (FR/RF) is a
    // valid inward-facing pair.
    if r1.strand == r2.strand {
        return PairRelation::Discordant;
    }
    if transcript_span(r1, r2) > PROPER_SPAN_MAX {
        return PairRelation::Discordant;
    }
    PairRelation::Proper
}

/// Span of the pair on the transcript: reference span from R1 start to R2 end
/// with the union of intronic RefSkip (N) intervals removed. RNA mates routinely bridge
/// kilobase introns; counting N toward the span would reject spliced proper
/// pairs that are adjacent on the transcript.
pub fn transcript_span(r1: &ReadAlignment, r2: &ReadAlignment) -> u32 {
    let left = r1.pos.min(r2.pos);
    let right = ref_end(r1).max(ref_end(r2));
    if right <= left {
        return 0;
    }
    let mut skips = Vec::new();
    for aln in [r1, r2] {
        let mut p = aln.pos;
        for op in &aln.cigar {
            match op {
                CigarOp::RefSkip(n) => {
                    let s = p.max(left);
                    let e = (p + *n).min(right);
                    if e > s {
                        skips.push((s, e));
                    }
                    p += *n;
                }
                CigarOp::Match(n) | CigarOp::Del(n) => p += *n,
                CigarOp::Ins(_) | CigarOp::SoftClip(_) => {}
            }
        }
    }
    skips.sort_unstable();
    let mut cursor = left;
    let mut n_total = 0u64;
    for (start, end) in skips {
        if end > cursor {
            n_total += u64::from(end - start.max(cursor));
            cursor = end;
        }
    }
    ((right - left) as u64).saturating_sub(n_total) as u32
}

/// One-past-the-end reference position: `pos` advanced by M/D/N ops.
pub fn ref_end(aln: &ReadAlignment) -> u32 {
    let mut end = aln.pos;
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => end += n,
            CigarOp::Ins(_) | CigarOp::SoftClip(_) => {}
        }
    }
    end
}

/// Signed template length from R1's viewpoint: positive when R1 is the left
/// mate, negative (magnitude = right-pair view) when R1 is on the right.
/// Magnitudes are transcript-aware (intronic N-skips excluded), consistent
/// with [`transcript_span`].
pub fn template_length(r1: &ReadAlignment, r2: &ReadAlignment) -> i32 {
    if r1.contig != r2.contig {
        return 0;
    }
    if r1.pos <= r2.pos {
        transcript_span(r1, r2) as i32
    } else {
        -(transcript_span(r1, r2) as i32)
    }
}

/// Integer-friendly insert statistics over |tlen| values (population
/// moments: `stdev = sqrt(max(0, E[x²] − m²))`).
#[derive(Clone, Copy, Debug, Default)]
pub struct InsertStats {
    /// Number of observations.
    pub count: u64,
    /// Σx.
    pub sum: f64,
    /// Σx².
    pub sum_sq: f64,
}

impl InsertStats {
    /// Record one |tlen| observation.
    pub fn push(&mut self, tlen: i32) {
        let x = tlen.unsigned_abs() as f64;
        self.count += 1;
        self.sum += x;
        self.sum_sq += x * x;
    }

    /// Mean of the observations (0 when empty).
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum / self.count as f64
        }
    }

    /// Population standard deviation (`sqrt(max(0, E[x²] − m²))`).
    pub fn stdev(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        let n = self.count as f64;
        let ex2 = self.sum_sq / n;
        let m = self.sum / n;
        (ex2 - m * m).max(0.0).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seed::Strand;

    fn aln(pos: u32, strand: Strand, cigar: Vec<CigarOp>) -> ReadAlignment {
        ReadAlignment {
            contig: 0,
            pos,
            strand,
            cigar,
            ..Default::default()
        }
    }

    #[test]
    fn opposite_strand_pairs_are_proper_in_either_order() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(1200, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Proper);
        assert_eq!(relate(&b, &a), PairRelation::Proper);
        assert_eq!(template_length(&a, &b), 350);
        assert_eq!(template_length(&b, &a), -350);
    }

    #[test]
    fn same_strand_pairs_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(1200, Strand::Plus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }

    #[test]
    fn intron_skips_do_not_count_toward_the_span() {
        // R1 bridges a 5 kb intron; the transcript-level span is 350 bp.
        let a = aln(
            1000,
            Strand::Plus,
            vec![
                CigarOp::Match(100),
                CigarOp::RefSkip(5000),
                CigarOp::Match(50),
            ],
        );
        let b = aln(6200, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(transcript_span(&a, &b), 350);
        assert_eq!(relate(&a, &b), PairRelation::Proper);
    }

    #[test]
    fn distant_mates_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(5000, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }

    #[test]
    fn different_contigs_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let mut b = aln(1200, Strand::Minus, vec![CigarOp::Match(150)]);
        b.contig = 1;
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }

    #[test]
    fn intron_union_shared_skip_is_subtracted_once() {
        for (strand1, strand2) in [(Strand::Plus, Strand::Minus), (Strand::Minus, Strand::Plus)] {
            let a = aln(
                100,
                strand1,
                vec![
                    CigarOp::Match(50),
                    CigarOp::RefSkip(2000),
                    CigarOp::Match(30),
                ],
            );
            let b = aln(
                120,
                strand2,
                vec![
                    CigarOp::Match(30),
                    CigarOp::RefSkip(2000),
                    CigarOp::Match(50),
                ],
            );
            assert_eq!(transcript_span(&a, &b), 100);
            assert_eq!(transcript_span(&b, &a), 100);
            assert_eq!(template_length(&a, &b), 100);
            assert_eq!(template_length(&b, &a), -100);
            assert_eq!(relate(&a, &b), PairRelation::Proper);
        }
    }

    #[test]
    fn intron_union_partial_nested_and_disjoint_intervals() {
        let cases = [
            // a skips[150,2150); b skips[170,2170): union2020,span100.
            (
                120,
                vec![
                    CigarOp::Match(50),
                    CigarOp::RefSkip(2000),
                    CigarOp::Match(50),
                ],
                100,
            ),
            // b's skip[170,2070) is contained in a's.
            (
                120,
                vec![
                    CigarOp::Match(50),
                    CigarOp::RefSkip(1900),
                    CigarOp::Match(50),
                ],
                80,
            ),
            // b skips[2200,2300), disjoint from a.
            (
                2170,
                vec![
                    CigarOp::Match(30),
                    CigarOp::RefSkip(100),
                    CigarOp::Match(50),
                ],
                150,
            ),
            // b skips[2150,2250), touching a's end.
            (
                2130,
                vec![
                    CigarOp::Match(20),
                    CigarOp::RefSkip(100),
                    CigarOp::Match(50),
                ],
                100,
            ),
            // Query-only ops do not expand the reference span; D does.
            (
                120,
                vec![
                    CigarOp::SoftClip(10),
                    CigarOp::Match(30),
                    CigarOp::RefSkip(2000),
                    CigarOp::Ins(3),
                    CigarOp::Match(40),
                    CigarOp::Del(10),
                ],
                100,
            ),
        ];
        let a = aln(
            100,
            Strand::Plus,
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(2000),
                CigarOp::Match(30),
            ],
        );
        for (pos, cigar, expected) in cases {
            let b = aln(pos, Strand::Minus, cigar);
            assert_eq!(transcript_span(&a, &b), expected, "{b:?}");
            assert_eq!(transcript_span(&b, &a), expected, "{b:?}");
        }
    }

    #[test]
    fn intron_union_overlap_cannot_make_excessive_span_proper() {
        let a = aln(
            100,
            Strand::Plus,
            vec![
                CigarOp::Match(50),
                CigarOp::RefSkip(2000),
                CigarOp::Match(30),
            ],
        );
        let b = aln(
            120,
            Strand::Minus,
            vec![
                CigarOp::Match(30),
                CigarOp::RefSkip(2000),
                CigarOp::Match(1150),
            ],
        );
        assert_eq!(transcript_span(&a, &b), 1200);
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
        let edge = aln(
            120,
            Strand::Minus,
            vec![
                CigarOp::Match(30),
                CigarOp::RefSkip(2000),
                CigarOp::Match(950),
            ],
        );
        assert_eq!(transcript_span(&a, &edge), 1000);
        assert_eq!(relate(&a, &edge), PairRelation::Proper);
    }

    #[test]
    fn intron_union_matches_independent_small_coordinate_mask() {
        // Enumerate small legal CIGAR geometries; a base mask supplies a
        // separate reference calculation without production merge helpers.
        for p in [0, 3, 8] {
            for q in [0, 4, 9] {
                for gap1 in [0, 2, 7] {
                    for gap2 in [0, 3, 8] {
                        let cigar = |gap| {
                            let mut ops = vec![CigarOp::Match(5)];
                            if gap != 0 {
                                ops.push(CigarOp::RefSkip(gap));
                            }
                            ops.push(CigarOp::Match(6));
                            ops
                        };
                        let a = aln(p, Strand::Plus, cigar(gap1));
                        let b = aln(q, Strand::Minus, cigar(gap2));
                        let end = (p + 11 + gap1).max(q + 11 + gap2);
                        let expected = (p.min(q)..end)
                            .filter(|&x| {
                                !(p + 5..p + 5 + gap1).contains(&x)
                                    && !(q + 5..q + 5 + gap2).contains(&x)
                            })
                            .count() as u32;
                        assert_eq!(transcript_span(&a, &b), expected);
                    }
                }
            }
        }
    }
}
