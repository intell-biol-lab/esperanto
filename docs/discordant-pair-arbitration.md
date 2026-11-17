# Experimental discordant full-body pair arbitration

`map --experimental-discordant-pair-arbitration` is default off. It applies
only to final paired output; discovery and older flow entry points preserve
their policy.

A full-body target can have MAPQ20 while its best independent placement is
discordant with an exact mate. The previous poor-placement trigger then skips
pair arbitration. A one-substitution true candidate can lose to an exact
remote homolog even when the pair supports the former. This option considers
already cached candidates under explicit fragment geometry, without adding
reference searches or changing sequence scores.

Both original placements must be full-length one-M paths, unrescued and
search-complete. The supporting mate must satisfy the existing pileup gate
and match the reference exactly. The target must still match its independently
scored cache. Editing-rescored pools are excluded. Candidate bodies must be
disjoint, inward-facing on the mate contig and within the existing proper-span
cap. Shared overlapping columns are not independent support.

The best compatible candidate must be a full-body ordinary genomic primary
path and cannot trail the old target by more than the existing UNIQUE_MARGIN.
The unchanged scorer retains within-fragment competitors; its computed MAPQ
must reach the existing pileup gate. New results are marked rescued and capped
at that gate. No MAPQ0 result is overridden upward. Conflicting proposals in
both directions are refused using the original pair, without order effects.

This is conditional reference-placement evidence. A mistaken supporting mate,
absent true source or missed candidate can still invalidate the conclusion.
The experiment must report per-read and per-pair gross losses, reference
equivalence/outside counts, true and foreign event observations at BQ13 and
known MAPQ20, fixed-catalogue no-edit false ALT, and DNA-confound controls.
Software tests do not prove confidence calibration, biological editing calls
or general superiority. Models, thresholds and the independent evaluator are
unchanged; integration, release and adoption require separate review.
