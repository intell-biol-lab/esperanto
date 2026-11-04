# Observable paired fragment endpoints — development proposal

This is a mapper development proposal, not scientific adoption or a release.
Its human replay, disabled ablation and independent source/family qualification
are still required. The model and independent evaluation contract are unchanged.

## Evidence before alignment

RNA paired reads default to an observable endpoint inference. DNA and single
reads do not use it. `map --no-fragment-endpoints` disables it; programmatic
callers can set `AlignConfig.fragment_endpoints=false`.

Require a common fragment length T >=31 and at least eight suffix bases in each
read. Each suffix must exactly match its corresponding builtin TruSeq or
Nextera adapter prefix, through the shorter of suffix and adapter length. Both
ends must support the same family. Bases after a complete adapter are not
interpreted. Enumerate the complete QC overlap domain (overlap >30 bases,
at most min(5,overlap/10) mismatches), using case-insensitive byte equality:
matching ambiguous bytes still count as matching when deciding whether an
offset competes. Reverse complementation preserves non-ACGT bytes unchanged,
exactly as QC does; alignment's unknown-to-N normalization must not hide them.

There must be exactly one viable offset. It must imply T=offset+R2_length,
have zero mismatches, and contain canonical bases only. The fragment must
contain at least three different bases and no base may exceed 80%. Missing,
short, inconsistent, erroneous, ambiguous or low-complexity evidence abstains.
Read names, reference placements, simulation labels and competitor output do
not participate.

## One effective alignment and durable raw records

Infer before chunk-local sequence interning. The key is the actual effective
sequence: supported T-prefix or unchanged raw bytes. Align it once through the
existing scorer, 0.66 coverage fraction and candidate/search budgets. Mate
arbitration uses those same effective bytes. There is no extra raw alignment.

For accepted primary and diagnostic candidates, restore the raw query frame
by adding raw_length-T soft-clipped bases at the sequenced 3' end: right on
plus, left on minus. Merge adjacent terminal S operations only. Preserve
reference coordinates, score, mismatch counts, junctions and MAPQ. Preserve
complete raw SEQ and QUAL with the existing strand convention.

This explicitly changes the conditional coverage denominator. An 80-base
supported interval requires 53 informative bases, not the raw 150-base read's
99. It is endpoint inference/internal preprocessing, not an isolated gain in
the extension algorithm. Compare raw input separately from equal-QC input.

Persist `XF:i:T` on mapped AND unmapped primary records. Valid XF endpoints
(31<=T<raw_length, BAM integer type) exclude the record from subsequent
raw-read collapsed rescue, deferred coverage rescue and stale-rescue merge
replacement. The durable raw BAM and mate slot determine exclusions on restart.
Excluded unresolved records remain in unmapped FASTQ even when unrelated reads
are rescued. The merge eviction re-reads that complete original pool rather
than exporting only eligible rescue entrants. Invalid XF is not valid evidence.

## Diagnostics and limitations

Directly affected rows use schema 2: raw read length/raw minimum remain present;
`effective_read_length`, `effective_minimum_primary_query_bases` and
`fragment_evidence` state the new eligibility explicitly. Candidate CIGARs
consume the raw sequence. Other rows retain schema 1. Two-pass discoveries
can indirectly affect schema-1 recipients; human comparisons must distinguish
direct endpoint changes from indirect junction-library effects.

Exact adapter-like biological observations can be indistinguishable from
read-through. A synthetic long biological molecule with identical read bytes
demonstrates 140 genuine body bases clipped across two 150-base reads when
T=80 is inferred. This is a real false-trimming limit, not a successful control
or promised origin discrimination. Extra independent library evidence is
needed to resolve it. The separate Myers neutral-terminal edit tie remains
unfixed; endpoint lifting does not establish terminal-edit retention.

Synthetic/CI success does not establish mapper superiority, calibrated MAPQ,
human qualification or adoption. Preserve every body/ALT loss and technical
gain in full development replay and future independent qualification.
