# Reference fragment decoding

Dev.38 changes `Contig::decode_into` and its existing append/slice callers.
It decodes packed bases over a contiguous range sharing one ambiguity decision.
Base access and all alignment search/scoring rules are unchanged.

Previously, each decoded base repeated `is_n`'s binary search over the same
interval table. `ambiguity_span` follows that same decision tree once. A
`pos < interval.start` branch sets an upper bound at `interval.start`; a
`pos >= interval.end` branch remains true for all larger positions; a match
remains valid until `interval.end`. Intersect those upper bounds with the
requested end. Before that bound every branch decision is unchanged, so the
same N/plain result applies to the whole span.

The implementation does not assume a sorted or merged table beyond the old
binary-search behavior. This matters because `Contig` fields are public and the
legacy index reader does not reject every malformed interval table. Independent
scalar-oracle tests include unsorted, overlapping and reversed intervals, all
small subranges, random arbitrary tables, appended prefixes and empty/reversed
requests. An operation-count regression checks that one 4,096-base plain span
no longer performs 8,192 N-table probes. The test-only counter is absent from
production.

Standalone decoding measurements establish only kernel cost. Adoption requires
unchanged independent scientific qualification and matched end-to-end runs.
