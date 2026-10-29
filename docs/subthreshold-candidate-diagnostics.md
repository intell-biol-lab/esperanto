# Executable candidate diagnostics contract

`esperanto map --candidate-diagnostics` writes
`OUT/candidate_diagnostics.jsonl`. The flag defaults to off. The sidecar is
evidence only: it does not supply primary/secondary BAM records, pileup inputs,
confidence calibration, or biological-origin labels.

Mapping policy, scoring, models and primary eligibility fraction (0.66) are
unchanged. For a 150 bp read, the minimum informative coverage is 99 bp:
positive legal 96/98 bp candidates remain visible but cannot become primaries.
An eligible fallback at another locus remains selectable, while a higher
scoring short competitor still contributes to existing ambiguity handling.

## Read rows

There is one JSON object per input read, in input order (mate 1 then mate 2 for
paired input). Identical sequences retain separate read identities. Pass1
discovery emits no rows.

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer `1` |
| `phase` | `pass2_before_collapsed_rescue` |
| `read_name` | Lossless UTF-8 FASTQ name, JSON escaped |
| `mate` | `1` / `2` for paired reads, `null` for single reads |
| `read_length` | Input query length |
| `minimum_primary_query_bases` | Ceiling of existing fraction times read length |
| `pass2_primary_placed` | Placement after pass2 mate arbitration |
| `independent` | Independent normalized candidate pool |
| `mate_attempt` | Last observed mate-window pool, or `null` if not attempted |

Collapsed rescue runs later in `map` and can change final BAM. Consumers must
join final BAM independently; pass2 status is not final mapping status.
Write/flush errors and invalid UTF-8 names fail the command. Partial output
from a failed command must not be used as complete evidence.

## Candidate pools

Every pool contains `stage`, `scored_candidates`, `excluded_nonpositive`,
`excluded_invalid_geometry`, and `candidates`. Exclusion categories are
disjoint: geometry is checked first, then score. The sum of exclusions and
retained candidates equals `scored_candidates`.

From dev.12, a limited mid-occurrence search adds optional `seed_budget`:
`eligible_anchors`, `retained_anchors` (at most 20,000), and `omitted_seeds`.
This limitation survives cache/mate arbitration and forces MAPQ 0. See the
[repeat-budget contract](repeat-budget-contract.md). Complete pools omit it.

Stages are `independent`, `mate_window_no_proposal` (empty attempt),
`mate_window_rejected` (after normalization, before the existing eligibility
check), and `mate_window_arbitration` (after the existing fragment reach
filter, before arbitration). Rejected attempts are retained. Pools are
retained heuristic proposals, not exhaustive alternative alignments. Earlier
deduplication, generator caps and filtering still apply.

Only positive common-score candidates with legal geometry are emitted:
nonzero M/I/D/N/S lengths, full query consumption, terminal soft clips, valid
reference/coordinate bounds, at least one M, and N bracketed by M on both
sides. Geometry exclusion is diagnostic validation only.

| Candidate field | Meaning |
| --- | --- |
| `rank_in_pool` | Zero-based position in original pool; exclusion gaps allowed |
| `reference_name`, `reference_start0` | Reference and zero-based start |
| `strand`, `cigar`, `origin` | Query alignment strand, reference-oriented CIGAR, existing generator label |
| `common_score`, `editing_score_mode` | Existing common score and `off` / `ag` / `tc` scoring mode |
| `ordinary_mismatches`, `edit_ag`, `edit_tc` | Existing common scorer counts |
| `informative_query_bases` | Existing informative query coverage |
| `eligibility` | `eligible` or `below_coverage` |

D and N do not contribute query coverage. Ambiguous columns follow the
existing scorer's treatment; M length alone does not establish eligibility.
Scores and counts are evidence, not calibrated probabilities or MAPQ.

## Executable checks

Mapper tests cover 96/98/99 bp in both orientations, ambiguous columns and
D/N coverage, geometry and score exclusions, different-locus fallback and
ambiguity, mate-window rejection/success/no proposal, deterministic ordering
across threads, duplicate sequence identities, pass1 suppression, unchanged
BAM/unmapped output, JSON escaping and IO/UTF-8 error propagation.

Full qualification also requires a new clean-build human replay, independent
scoring/audit, and diagnostic-on/off comparison of all final primary records.
Synthetic tests and CI alone do not qualify a candidate for scientific use.
