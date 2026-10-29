# Candidate coverage, deduplication and fallback contract

This development contract applies after candidate generation from genomic,
transcriptome and mate-window paths. It separates the score of an explanation
from its eligibility to become a primary read placement. Its mechanism tests
are synthetic and do not establish general mapping accuracy.

## Coverage eligibility

For a read of length `L > 0`, let `Q` be informative query bases in M columns
plus query bases in I columns. An M column containing N on either side does
not contribute, following the existing scorer. D and N CIGAR operations
consume reference, so neither contributes to `Q`. S contributes nothing.

Coverage is sufficient exactly when `Q >= 0.66 * L`, using the existing
`REFUSE_MIN_ALIGNED_FRAC`. For a 150-base read, 98 fails and 99 passes;
plus and minus orientations use the same rule. Coverage does not mean exact
matches: informative mismatches count toward coverage and affect the score.

`StrictScore.aligned` historically includes deleted reference columns and
feeds edit-density consumers. This change leaves that bookkeeping and all
scoring constants intact. The coverage predicate subtracts CIGAR D length.
It does not use reference span, intron length, producer score or chain score.

## Deduplication and editing evidence

Candidates use the existing deterministic decision order: score descending,
ordinary mismatches ascending, origin rank, contig, position. The existing
same-locus footprint definition (same contig/strand, reference M overlap at
least 0.5 of the smaller footprint) is unchanged.

Deduplication may remove a candidate only in favor of a better-ranked
candidate in the **same coverage class**. A short high-scoring explanation
therefore cannot remove a lower-scoring full-coverage explanation.
Nonpositive full-coverage candidates are retained at this stage because a
single-direction editing signature may invoke the existing reduced-alphabet
rescoring and make them positive. Hyper rescoring uses the same coverage-class
dedup after extension and score updates; eligibility is never cached as a
stale flag. Direction thresholds and scoring remain unchanged.

## Primary selection, fallback and ambiguity

A candidate can become primary only with sufficient coverage and a strictly
positive score in the active scoring domain. Scan the ordered candidates to
the first such candidate. If none qualifies, return unmapped so downstream
rescue can act. A full-length negative or zero score is not a fallback.
Do not force an old wrong alignment back into the output to recover a metric.

The runner-up is the highest positive-scoring **different-locus** explanation
among all retained candidates, including those with insufficient coverage.
A stronger short explanation at another locus is unresolved ambiguity:
its score can exceed the selected fallback and produce MAPQ 0 under the
existing formula. Same-locus short explanations do not create a second locus.
No MAPQ formula or ceiling changes are introduced.

`ReadAlignment.n_cands` counts retained **primary-eligible** locus
representatives; coverage-ineligible evidence can affect `second_score`
without increasing that count. `genomic_any` retains its existing meaning
that the genomic generator produced evidence, rather than being a new
qualification flag.

## Mate-window path

Mate re-arbitration uses the same coverage classes and primary predicate.
Check the first eligible candidate when testing whether the mate-window
locus won; an empty eligible pool returns `None` safely. Pair-distance
pruning and the rescued MAPQ ceiling remain unchanged. Cached independent
candidates preserve both coverage classes for this path.

## Executable checks and scope

Eight `coverage_*` tests in `crates/map/src/align.rs` cover:

- A high-scoring same-locus short hit cannot delete the full fallback.
- A disjoint short score leader permits fallback while preserving MAPQ ambiguity.
- Nonpositive full candidates cannot become primary.
- The 98/99 boundary on both strands and insertion contribution.
- Deletion/skip exclusion and ambiguous M columns.
- Full negative editing evidence survives dedup and can become eligible.
- Mate-window rejection with no eligible candidate, and a valid mate positive.
- Dedup and fallback remain independent of generator order.

Five mechanisms fail against dev.3 behavior before this change. Full tests,
strict human interface qualification and scientific per-read/event review
are separate gates. Search remains heuristic and the footprint relation is
not a proof of sequence equivalence. This contract does not relax the
coverage fraction, retrain a model or modify the independent evaluator.

## Intron fast-path candidate independence (dev.5)

An intron-chain fast proposal searches the complete anchor pool. Its exon pair
can differ from the chain that invoked it. Do not replace its junction with
that chain's boundary or change only its N length: that does not move the
CIGAR donor or query breakpoint. Before retaining a fast proposal, reconstruct
each skipped reference interval from position and CIGAR and require equality
with the declared junction intervals, including contig and multiplicity.

Whether a fast proposal succeeds or is rejected, continue the independent
generic spliced path. Both proposals enter existing common rescoring and
coverage-class dedup. This restores candidate availability without forcing a
particular locus or lowering eligibility. Synthetic tests cover both strands,
CIGAR/query-reference consumption, missing and inconsistent junction metadata,
and an equal-locus-pair editing trigger that must retain MAPQ 0 ambiguity.
This does not guarantee recovery of an exon with no chain anchors or calibrate
MAPQ against a population of reads.
