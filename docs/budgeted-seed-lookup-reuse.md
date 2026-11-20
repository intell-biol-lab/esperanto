# Budgeted seed lookup reuse

The budgeted ordinary anchor collector previously queried each read minimizer
for eligibility, then queried eligible minimizers again for position emission.
The complete-budget path repeated queries even for absent/above-cap minimizers.

Retain the immutable `IndexHit` from the first query instead. This does not
copy posting lists or create cross-read caches. Each input minimizer is still
queried separately, including repeated k-mers at different query positions.
The complete path emits eligible lists in original input order. The overflow
path keeps the existing `(count, pos, kmer, strand)` ordering and whole-list
admission, with unchanged omitted-seed accounting. No budgets, candidates,
reference layout, scoring, tie policy, MAPQ or model parameters change.

The per-call eligible-seed scratch entry grows from `(index, count)` to
`(index, IndexHit)`; this trades a larger temporary descriptor for eliminating
binary searches. Posting storage remains borrowed. Overflow sorting may pay
for the larger descriptor, so microbenchmark both complete and overflow paths.

Tests include a lookup-count regression (red on the parent behavior), inherited
budget/order/strand tests and a frozen-parent differential oracle with 22,528
budget/cap/permutation cases. These establish collector equivalence on tested
inputs, not universal mapper qualification or an end-to-end speedup.

Research evidence lives in the independent ESP research workspace under
`experiments/esp-seed-lookup-reuse-20261120-attempt01/`. Production binaries and
software/model adoption locks are not changed by this candidate.
