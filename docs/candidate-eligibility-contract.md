# Candidate coverage, deduplication and fallback contract

This development contract applies after candidate generation from genomic,
transcriptome and mate-window paths. It separates the score of an explanation
from its eligibility to become a primary read placement. Its mechanism tests
are synthetic and do not establish general mapping accuracy.

## Independent full-anchor intron proposals

For RNA, a best chain with at least two segments may invoke the existing
two-fragment intron constructor on the complete same-contig/strand anchor
pool. The chosen pair need not belong to that chain. Three or more segments
do not disqualify a verified pair. Preserve the constructor's verification
and geometry checks, generic spliced proposals, continuous segment competitors,
common scoring and primary coverage gate.

Executable regression:
`intron_pool_proposal_survives_three_segment_best_chain` covers both read
alignment strands, either terminal clip, pre-existing full-pool support,
proposal generation before dedup, continuous competition and arbitration.

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
ordinary mismatches ascending, informative query coverage descending, origin
rank, contig, position. The coverage tie preference was added in dev.7. The existing
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

## Fast proposal reconstruction (dev.6)

The fast constructor retains its existing fragment-pair selection and dual
mismatch split-verification gates. Reconstruct from that selected pair: refine
its own tight read boundary with the existing read-aware routine (or retain
the pair's diagonal boundary when the routine cannot apply), concatenate
flanking exon windows, extend against that pseudo-reference and stitch using
the very same junction. The emitted score is the real extension score, rather
than a fabricated full-read match total. Retain legal outer clips and indels.
Do not rebuild from an unrelated invoking chain. The generic spliced path
still runs independently. A synthetic asymmetric 60/79-base splice with an
11-base technical clip fails construction consistency on dev.5 and passes
on both read strands and both genomic motif orientations after reconstruction.

## Neutral endpoints and boundary competition (dev.7)

After extension traceback/end trimming, prefer the longest available contiguous
diagonal outer flank whose net substitution score is exactly zero. A positive
alignment is required; ambiguous bases on either side stop the scan. This is
a tie convention, not a clipping penalty or editing exemption. Existing
substitution/gap scores, retained internal indels and total score stay unchanged.
The same shared finish routine is used by scalar, AVX2 and banded kernels.
This does not claim a globally longest optimal gapped alignment.

On equal common scores and ordinary mismatch counts, informative query coverage
breaks candidate ties before origin and coordinates. D and N cannot inflate
this quantity. Both generator orders must keep a same-locus terminal edit when
the clipped explanation only ties its score. Different-locus competition and
the primary coverage/positive-score rules are unchanged.

Read-aware splice refinement must not equate a mate's alignment direction with
the RNA direction. Search each existing boundary shift under both transcript
orientations with the same existing flank/motif/library scoring. Preserve the
selected transcript's junction metadata and the read's own alignment strand
separately. Known and novel canonical junctions are available on both mates.
This matches the existing any-strand common junction scorer; it does not add
an annotation-only override or force a truth-derived boundary.

Five synthetic tests cover neutral Ag/Tc flanks at both ends, negative/N
controls, input-order invariant same-locus ties, annotated/novel 1-2 bp
boundary competition on both motif/read orientations, and an end-to-end
spliced terminal-edit competition with both mates. Four fail on dev.6.

## Ordinary-mismatch barrier for neutral outer flanks (dev.8)

The dev.7 endpoint convention omitted the ordinary-mismatch part of common
arbitration. A zero-net-score extension with additional ordinary mismatches
must not beat clipping just by covering more query bases. Stop an outer-flank
scan at an ordinary mismatch or an ambiguous base; never cross that barrier
to find a more distant zero-net-score endpoint. Keep a valid neutral edit
block reached before the barrier. Exact matches and directional A>G/T>C
columns use the shared `exempt_dir` classification, independently of which
extension pass currently exempts one editing direction.

This changes only the shared finish routine's tie normalization. Positive
score extensions, substitution/gap constants, coverage eligibility, splice
refinement, model and evaluator remain unchanged. No adapter identifier,
biological-origin label, human coordinate or new length/penalty threshold
enters the mapper.

Synthetic non-N adapter, poly(A) and deterministic random-tail competitions
cover left/right ends, Off/Ag/Tc scoring, and scalar/packed dispatch. Separate
controls preserve terminal editing and positive-score ordinary alignments.
Identical sequence/reference inputs with different body/technical labels
remain indistinguishable: sequence-only tie handling cannot guarantee perfect
technical-tail specificity. Human replay and independent held-out controls
are separate gates.
