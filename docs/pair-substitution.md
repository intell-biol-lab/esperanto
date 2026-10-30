# Pair-supported sparse-substitution endpoint competition

Dev.23 is a development proposal, not a calibrated mapping or biological
editing claim. It runs only on paired RNA with existing exact read-through
endpoint/overlap evidence, in the final output pass before raw CIGAR lifting.
`map --no-pair-substitution` independently disables it. Older flow entry points
delegate to `map_stage_with_pair_substitution_options` with the RNA default.
DNA, single-end and transformed rescue do not use this pair transaction.

Both original mapped mates must be nonrescued, not search-limited, MAPQ20..254,
same contig/opposite strand and proper by union-intron span <=1000. Their
original paths must be M/S only with at least one effective terminal clip.
All aligned nonreference columns must belong to one A>G or T>C class.
At least two distinct aligned reference positions must carry the same
nonreference base in both mates at Phred >=20. Unknown or ordinary/mixed
substitutions refuse the witness. Phred and sequence are oriented from the raw
effective query for each mate; a mate is not an independent biological molecule.

That witness defines one read-oriented scoring domain per mate. All cached
proposals and the original alignment are retained alongside bounded unspliced
extensions. Every candidate is rescored in that domain, complementing A>G/T>C
for a candidate on the opposite strand. Existing coverage-aware footprint
deduplication and ordering apply; this is not complete genome search.
The separately rejected exact-path arbitration policy stays off by default.

The winning ordinary coverage-qualified proposal must strictly improve the
original score in the conditioned domain, retain the same contig/strand and
every old aligned query-to-reference column, add coverage and contain only
matches or supported substitutions in newly added columns. M/S geometry is
required for acceptance. A different winning competitor refuses the transaction;
there is no fallback to a favorable lower-ranked extension. Neutral score ties
preserve the original alignment. Both decisions use the original pair snapshot.

Accepted records carry MAPQ0 and `RE:Z:pair_substitution`; the real conditioned
scores and reference mismatch counts are preserved. They do not enter pass1
junction discovery or raw collapsed re-rescue. Candidate diagnostics report
the reference-oriented witness positions, direction, original conditioned
score, full scored pool and acceptance decision without altering output.

Shared substitutions establish a mapping scoring condition, not editing
biology, exhaustive uniqueness or probability calibration. Missing true
sources and adapter-like biological ends can produce identical observations.
No claim to distinguish such worlds follows. Full replay must account for
gross body/ALT losses, technical gains, wrong placements and MAPQ gate costs;
source/region/family holdout and fair broader competitor evaluation remain
necessary to establish mapper advantage.
