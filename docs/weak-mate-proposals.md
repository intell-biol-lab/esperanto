# Experimental weak-mate candidate proposals

Prospective development contract, 2026-10-30; disabled by default.
The independent source-isolated development panel exposed a candidate omission:
one read's retained candidates omit its ancestral repeat copy, while its mate
has an exact ordinary M/S path there but MAPQ0 after an incomplete seed search.
This is a candidate generation issue; MAPQ0 is not evidence of a unique locus.

An opt-in experiment may use that mate solely to propose one bounded ordinary
window alignment. Both original reads must already be mapped and unrescued.
The target must be below the pileup gate; the proposer must have MAPQ0, disclosed
seed-search incompleteness, an exact primary-eligible M/S path and positive score.
Validate proposer reference/query geometry and bases; do not trust score alone.

The additional target must be exact, primary-eligible, gap-free ordinary M/S,
proper relative to the proposer, and strictly exceed the original target and
every retained distinct path after common ordinary rescoring. Keep every cached
competitor, including other chromosomes, coverage-ineligible paths and explicit
original geometry. Never prune competitors by pair consistency or footprint.
Any equal-or-better distinct path refuses replacement. No uniqueness claim is
made: all accepted replacements are rescued and forced to MAPQ0, with incomplete
search disclosed. Candidate diagnostics use stage `weak_mate_proposal`.

Both directions use the same original pair snapshot; never chain proposals.
The experiment runs only in pass2; it cannot alter the discovery library.
Existing confident-mate rescue and unmapped short-body policy stay separate.
An existing editing-domain candidate pool is refused by this ordinary proposal.
This does not alter editing domains, models, thresholds or evaluator semantics.

Synthetic controls cover weak versus confident proposers, exact and mismatched
reference paths, strand/geometry, strict improvement, equal/stronger remote
competitors, duplicate geometry, ablation and immutable caches. Full human
qualification compares old accepted development inputs and fresh development
inputs, raw and common-QC, both seeds. Reject adoption on any body/ALT loss or
new technical alignment at an original query offset. Test outcomes remain sealed.
Default-off parity with frozen dev.26 is an additional gate. This experiment
does not authorize merge, release or adoption-lock changes.
