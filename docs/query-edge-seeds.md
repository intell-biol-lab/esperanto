# Bounded query-edge seed completion

Development candidate: `1.1.1-dev.16`. Not a release or evidence of general
mapper superiority.

## Mechanism

Reference and read minimizers still use the original k/w sketch. A clipped
query can lose a minimizer selected with longer context even when that k-mer
is wholly inside the retained sequence. In RNA mode with editing extension
off, after both ordinary chain rounds fail to produce a chain, complete the
two query edges with valid canonical k-mers from the first and last `w - 1`
starts. For k15/w5 this adds at most eight seeds. No outside-query bases,
unknown-containing k-mers, source labels or expected alignment coordinates
are consulted.

Deduplicate the union by query start/canonical k-mer, sort deterministically,
and replace the failed mid-occurrence round using that union. The existing
empty-sketch early refusal is unchanged. A nonempty successful first or mid
round bypasses this fallback, as do DNA and editing-variant modes.

## Unchanged limits and semantics

- Reference index, sketching, occurrence cap 5000 and whole-list anchor budget
  20000; no coordinate-prefix sampling.
- Eight candidate chains and two shared DP calls; chain, coverage and scoring
  thresholds unchanged.
- Limited-search evidence remains attached through primary and mate
  arbitration, with MAPQ zero.
- Fragment endpoint inference and its raw SEQ/QUAL/terminal-S contract,
  collapsed rescue exclusion, classifier model and independent evaluators.

The fallback adds a bounded seeding/chaining attempt, so unchanged DP/anchor
limits are not a claim of unchanged total execution time.

## Executable evidence

`seed.rs` tests canonical encoding/orientation, both edges, overlapping ranges,
deduplication, idempotence, invalid domains/bases and stale lengths.
`align/prefix_tests.rs` constructs a wholly synthetic repeat/background/target
reference with the native frequency policy: raw context maps while its exact
prefix has no ordinary chain. Both strands recover without changing the
reference sketch. Separate controls exercise successful-round and other-mode
bypass, negative terminal score, and real overflowing repeat pools with MAPQ0.
`pipeline/fragment_integration.rs` checks actual paired BAM placement, terminal
soft clips, raw sequence/quality and the unique integer XF endpoint tag.

These controls establish mechanism behavior only. Independent full development
replay, every gross adverse body/ALT/technical column, ordinary controls,
competitor complements, source/family holdout and resource/confidence validation
remain separate scientific gates. This patch does not repair unrelated
terminal-edit scoring or biological/technical origin ambiguity.
