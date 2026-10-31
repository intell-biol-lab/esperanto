# Continuous extension windows for sparse anchors

The prior continuous window padded both reference ends by `read_length/2+30`.
When anchors covered only one end of a read, this could exclude the other end's
true reference bases before either verifier saw them. A synthetic 150-base
exact match with query anchor span `[2,42)` produced `147M3S`.

Dev.14 shares one window function between the short-read Myers verifier and
banded DP. Head padding is the greater of the prior flank and query span start;
tail padding is the greater of that flank and uncovered query suffix. Reference
bounds, the existing `ref_start+read_length+flank` right cap, and maximum width
`read_length+2*flank` bound the interval. Additions saturate.

For a collinear chain whose full-read diagonal lies inside the reference,
the interval includes the full diagonal. If both query overhangs fit within
the old flank, the interval is exactly the old interval. This does not promise
exhaustive indel search. Fast-lane eligibility, scores, candidate scheduling,
coverage and the two-DP budget are unchanged.

Six synthetic tests cover both strands, anchors near either end, short and
long reads, contig edges, N/non-N technical tails against concrete mismatching
reference flanks, DP terminal edits,
zero DP budget, displaced anchors, overflow boundaries, full-diagonal inclusion
and exact preservation of unaffected windows.

## Separate endpoint limitation

The new controls exposed an existing short-read Myers endpoint tie rule:
zero-score terminal editing substitutions can still be clipped by its Kadane
scan, even when the window contains them. DP has separate neutral-flank
retention logic. The terminal-edit test explicitly forces short reads to DP
with an internal ambiguous query base outside its exact anchors, and also
tests unambiguous long reads. It does not claim Myers endpoint parity.

The first mixed-path failure and its source are retained in the research
dev.14 evidence directory for a separate mechanism iteration. Positive-score
accidental adapter matches, MAPQ calibration, independently held-out performance
and resource comparison also remain open. A full development replay must
measure changes in body, editing observations and technical columns before
any adoption decision.
