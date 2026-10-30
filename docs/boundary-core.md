# Shared alignment columns at an uncertain homopolymer boundary

Dev.26 adds a final paired-RNA fallback, independently disabled with
`map --no-boundary-core` or `AlignConfig.boundary_core=false`. The existing
short-body enable switch remains required. DNA/single-end/coordinate-only mate
APIs cannot use it. Original short-body qualification runs first: any accepted
or ambiguously qualified original decision is preserved.

An exact contiguous below-coverage M/S candidate may have one pure canonical
terminal S>=30 whose base continues into adjacent M columns. Those matching
run columns cannot be assigned to biological versus technical sequence from
these observations alone. The old exact-boundary helper continues to refuse
certifying such a boundary. Instead, derive the alignment columns common to
all boundary assignments by leaving the entire maximal adjacent run unaligned.
Advance reference start for a prefix run; shorten only the final M for a suffix
run. Do not select an interior boundary, infer technical ancestry, modify raw
sequence/quality, or change an already mapped read.

The core must independently pass the existing exact/complex-body >=31 rules,
nonrescued MAPQ20..254 proper mate, global or mate-conditioned index support,
whole-list budgets and exact-path competition. Indels, splices, mismatching
source columns, empty/noncomplex cores and multiple supported core paths refuse.
All accepted projections remain MAPQ0 / mate_short_body and do not enter the
ordinary MAPQ20 pileup. Real biological run bases can also remain clipped:
the projection guarantees no biological/technical segmentation or completeness.
A wrong mate, missing source, omitted candidate or incomplete sketch remains
an unresolved source-location risk.

Diagnostics retain the original candidates and `boundary_core_proposals`:
source/core ranks, `uncertainty_query_start`, `uncertainty_bases`, `run_base`
and `assignment=unknown_biological_or_technical`. Coordinates use the declared
raw alignment query orientation, including minus-strand effective-fragment
lifting. Exact core paths deduplicate only after support qualification.
If several original extensions project onto the same selected core, preserve
all source records and emit the union of their uncertainty intervals.

Final BAM tag `XB:Z:query_start=Q;length=N;assignment=unknown` records that union,
using a zero-based start in SAM sequence/query orientation. These columns are
soft-clipped. `XB` never labels the rest of a soft clip as technical sequence.
The existing `RE:Z:mate_short_body`, evidence tags and raw SEQ/QUAL remain.

Native controls cover identical observed worlds with all possible run-boundary
assignments, four run bases, both strands/ends, runs longer than development
examples, source mismatch/I/D/N refusal, nearby omitted repeat copies, merged
core provenance, uncertainty/seed frame lifting, mate order, threads and one/two
passes. Fully biological already mapped terminal runs remain mapped unchanged.
The disabled arm must reproduce the exact parent mapping/scoring content in full
human replay. Software checks, partial read recovery and this invariant-column
argument do not prove calibrated confidence or general mapper superiority.
