# Mate-conditioned primary-index support

Dev.25 adds a bounded support route inside the existing paired short-body
qualification. The global-index support route remains unchanged. Only otherwise
unplaced paired RNA reads can enter; `--no-mate-short-body` disables both routes.
Use the exact dev.24 binary as the conditional-mechanism comparison arm rather
than interpreting that broader switch as a conditional-only ablation.

The new route requires the same exact complex body, positive score, below-0.66
ordinary coverage, a single canonical pure terminal S>=30 with a different
adjacent body base, and the existing nonrescued MAPQ20..254 proper mate. Only
unspliced M/S candidates use conditional support; global support still allows N.
The ordinary mapped path, model, thresholds, scoring and evaluator do not change.

Scan complete native body minimizer posting lists in ascending
(count, query position, canonical key) order. Each list has at most 5,000
postings and each candidate attempt scans at most 20,000 total postings, reusing
ordinary mid-occurrence limits. Skip oversized lists in full. For each posting,
derive the contiguous body's hypothetical start from the seed's body offset.
Require supplied-reference bounds, matching canonical orientation, mate contig,
opposite read strand and the existing <=1000 union-intron proper-pair span.
A list supports only when exactly one compatible start remains and it equals
the candidate start. Two disjoint query seed intervals must support the path.
Multiple qualifying exact paths refuse, including duplicates resolved only
after qualification.

This is mate-conditioned sketch support, not global uniqueness, exhaustive
genomic search, a likelihood or a calibrated probability. A wrong mate, missing
source or unrepresented competitor can remain misleading. Nearby indexed copies
defeat support even when full candidate generation omitted them. Matching
biological homopolymer boundaries remain refused. Accepted placements retain
MAPQ0, `RE:Z:mate_short_body`, and exclusion from the ordinary MAPQ20 pileup.
One rescued mate cannot bootstrap another.

Diagnostics add `support_scope` (`global_index` or `mate_conditioned_index`).
The old `disjoint_unique_indexed_seeds` contains only global support;
`disjoint_mate_conditioned_indexed_seeds` records conditional support separately.
Both use the declared raw alignment query frame, with effective-fragment lifting.
`conditional_sketch_attempts` records candidate rank, scanned lists/postings,
omitted and over-cap lists, and retained disjoint seed count, including failures.
Costs are per candidate and include repeated lists if repeated candidates are
attempted. Omitted lists were not searched; MAPQ0 never implies completeness.

Native tests cover remote versus nearby repeats, omitted nearby candidate paths,
both orientations and clip ends, raw-frame lifting, whole-list budgets and caps,
mate order and one/two-pass pipeline behavior. Existing exact-boundary,
confidence/span, indel/splice and biological-run controls still apply.
Full human replay and independent source/region/family generalization are
separate scientific gates. This proposal does not establish mapper superiority.
