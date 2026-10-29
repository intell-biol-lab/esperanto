# Bounded repeat candidate search

This development proposal replaces the mid-occurrence round's all-or-nothing
20,000-anchor abort. It does not change the 5,000 per-seed occurrence cap, the
round trigger (best round-one query-span coverage below 0.6), extension budget,
strict scoring, 0.66 primary coverage, or model.

Preflight the hit counts for every passing minimizer. When their total fits
20,000, preserve the existing collector's exact output and order. On overflow,
sort seeds by `(hit count, query position, canonical kmer, strand)`. Admit an
entire seed's hit list only if it fits the remaining budget. Omit the entire
list otherwise, continuing with later seeds. Never sample a coordinate prefix.
The full-anchor intron constructor uses the same bounded collection when
this search was limited. This contract bounds the mid-occurrence round and
that constructor; it is not a new global bound on the inherited first round,
editing-variant search, or total mapping runtime.

Record `eligible_anchors`, `retained_anchors` and `omitted_seeds` for a limited search. Propagate
that status even through an empty independent candidate cache, L1 winners,
mate-window fallback and later confidence overrides. Every resulting placed
read has MAPQ 0. It cannot establish genomic uniqueness even if only one
retained candidate remains. Retained repeat placements may still be wrong.
Unplaced limited reads retain the limitation in raw BAM EK evidence, without
changing FASTQ names. Collapsed rescue joins this evidence by exact template
name and mate slot, combines it with the transformed search's limitation, and
retains it in the final BAM, including interrupted merge/eviction reruns.

The existing EK fields keep their order. Limited records use
`mapq_src=seed_budget` and append `seed_search=limited`. JSONL candidate pools
carry an optional `seed_budget` object with the three counts above. The sidecar
remains before collapsed rescue; join the final BAM to determine final status.
Complete-search records keep the existing format. Consumers must tolerate the
documented optional fields.

Executable controls in `index::budget_tests` cover complete-order equivalence,
whole-list inclusion, budget boundaries, rare-seed preference and deterministic
overflow ordering. The `align` tests cover 700 identical loci on both alignment
strands with either technical terminal clip, single-path confidence, later
overrides and an empty-cache mate fallback. These tests must pass with the
proposal and the repeat test must fail with the old whole-round abort.
Flow controls cover rescue inheritance, mate separation and durable reruns;
both BAM encoder paths must preserve unplaced evidence.

Human development replay is a separate gate. Reject adoption on any added
per-read body or ALT loss, increased technical alignment, budget violation,
or limited-search MAPQ above zero. Passing software tests does not establish
scientific qualification or replace an independently frozen holdout.
