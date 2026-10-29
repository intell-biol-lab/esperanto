# Changelog

## [Unreleased]

### Organization development

- In `1.1.1-dev.18`, require additional ordinary repeat-context proposals to
  beat the strict local score against constant A/C/G/T backgrounds. Retain old
  candidate families and conservative MAPQ0; declare background score/rejection
  counters and independent `map --no-context-background` ablation. Generated
  both-strand and paired failures turn green while full low-complexity target
  controls remain. Independent human qualification and adoption are separate.

- In `1.1.1-dev.17`, propose bounded high-frequency ordinary context candidates
  after original chain failure. Preserve existing candidates/DP,verify legal
  local M/S hypotheses and propagate typed additional posting-budget evidence
  through the existing conservative MAPQ0 policy. Add independent
  `map --no-repeat-context` ablation; converted rescue remains unchanged.
  Generated AG/transversion both-strand failures reproduce on the parent.
  This remains unadopted development pending full scientific qualification,
  confidence and independent heldout evaluation.

- In `1.1.1-dev.15`, propose paired read-through endpoint inference from exact adapters and a
  unique QC-domain overlap before sequence interning/alignment. Use explicitly
  supported effective intervals, preserve raw SEQ/QUAL/CIGAR consumption, and
  persist XF on placed/unplaced reads. Schema-2 diagnostics declare conditional
  coverage. Exclude XF-bound reads from later raw collapsed rescue and preserve
  excluded unresolved reads during unrelated rescue/eviction and restart.
  This changes endpoint preprocessing/eligibility; it is not human-qualified,
  adopted or released. See the paired-fragment endpoint contract and its
  demonstrated biological adapter-mimic limitation.

- In `1.1.1-dev.9`, add opt-in `map --candidate-diagnostics` JSONL evidence
  for positive legal independent and mate-window candidates, including those
  below primary coverage. Preserve input/mate identities and pass2 ordering;
  label the phase before collapsed rescue. Mapping, MAPQ and the 0.66
  coverage policy are unchanged. See the executable diagnostics contract.

- Adopt `intell-biol-lab/esperanto` as the continuing development repository,
  preserving upstream history, attribution and release tags.
- Identify development builds as `1.1.1-dev.1`, with repository, source commit
  and dirty-state information in `--version`. The software version change
  intentionally prevents silent resume of runs stamped with upstream 1.1.0.
- Pin Rust 1.94.0 and establish CI for formatting, locked workspace tests,
  portable CLI release builds and bundled model integrity.
- Disable self-update for development builds and replace the inherited binary
  installer with source installation guidance. Future stable update endpoints
  name the organization; stable release hardening is still required.
- Record source/model adoption identities, contribution / release contracts
  and the unresolved internal-soft-clip qualification failure. No mapper
  repair or model change is included.

### Changed

- Updated the human model bundle. The classifier is retrained on an expanded, donor-partitioned panel (64 training, 16 calibration, and 20 held-out evaluation donors), and the decision threshold is recalibrated on the calibration donors. The 27-feature contract, input requirements, and output formats are unchanged.

### Fixed

- In `1.1.1-dev.16`, complete query-edge canonical seeds only after both
  ordinary RNA chain rounds fail. Replace the failed bounded mid round with
  the original/edge seed union, without changing the reference sketch,
  occurrence/anchor limits, extension budget, scoring or successful rounds.
  Retain limited-search MAPQ0. See the query-edge seed contract; this is not
  scientific adoption or a general mapper superiority claim.

- Remove two unjustified full-matrix `read_length<=reference_length` debug
  assertions. Short splice pseudo-references are valid rectangular local
  alignments. Keep release recurrence, widths, dispatch and scores unchanged;
  expand scalar/AVX2 rectangular parity controls. See the local-window contract.

- In `1.1.1-dev.14`, size both continuous extension windows using uncovered
  query ends as well as the existing flank. Sparse anchors no longer truncate
  a physically available full-read diagonal. Share geometry between Myers
  and DP while preserving the prior maximum width, right cap and DP budget.
  See the sparse-anchor window contract; independent qualification is required.

- In `1.1.1-dev.13`, preserve the nonempty complete first round's extension
  priority when repeat-seed search is limited. Other rounds still fill remaining
  slots and compete by common sequence score. Keep the eight-chain shortlist,
  shared two-DP budget, coverage, scoring, model and limited-search MAPQ 0.
  See the limited-search extension-priority contract. Independent qualification
  and per-read regression review remain required.

- In `1.1.1-dev.12`, replace the mid-occurrence round's whole-round abort
  with deterministic whole-seed-list selection within the same 20,000-anchor
  budget. Preserve complete-input order, cap limited-search placements at
  MAPQ 0 through cache/mate/override paths, and expose the limitation in EK
  and optional diagnostic evidence. See the repeat-budget contract. This
  remains a development proposal subject to per-read regression review.

- In `1.1.1-dev.11`, allow the existing full-anchor intron-pair proposal
  when the best chain has at least two segments. Extra repeat-derived
  segments no longer suppress an independently supported two-exon path.
  Keep fragment verification, continuous competitors, common rescoring,
  coverage, anchor/extension budgets and model unchanged. Synthetic cases
  cover both alignment strands and either terminal clip.

- In `1.1.1-dev.10`, keep hyperediting re-arbitration in one original-read
  conversion domain: complement A>G/T>C when a candidate's alignment strand
  differs from the signature-producing candidate. Apply that direction to
  genomic re-extension and common rescoring. Preserve signature thresholds,
  score constants, primary coverage and model weights. Four synthetic tests
  cover both strands/directions, L1/genomic paths and a same-strand negative
  control. See the strand-aware editing-domain contract; development replays
  do not establish independent qualification or calibrated MAPQ.

- In `1.1.1-dev.8`, apply the existing score/ordinary-mismatch decision order
  to score-neutral outer-flank normalization. Stop before ordinary mismatches
  instead of extending through them to a later score tie. Preserve A>G/T>C
  terminal evidence using the existing substitution classification. This
  prevents neutral non-N technical tails from buying coverage at the cost
  of extra ordinary mismatches. Positive-score tails and identical-sequence
  origin ambiguity remain separate qualification limits.

- In `1.1.1-dev.7`, retain score-neutral informative diagonal outer flanks
  after local extension, preserving terminal editing evidence on score ties
  while keeping negative and ambiguous clips. Same-score, same-mismatch
  candidate ties prefer informative query coverage. Read-aware junction
  refinement considers both transcript orientations, matching unstranded
  common arbitration instead of equating mate direction with RNA direction.
  Score constants, density/coverage thresholds and model weights are unchanged.

- In `1.1.1-dev.6`, reconstruct fast intron candidates by extending a
  pseudo-reference from their own fragment pair, then stitch at the same
  refined junction. This replaces the midpoint M-N-M assembly that could
  disagree with refined genomic boundaries and drops neither outer clips nor
  indels. Existing split-verification thresholds remain unchanged. The dev.5
  seam rejection alone lost supported fast candidates in human replay;
  that failed candidate and its results remain historical evidence.

- In `1.1.1-dev.5`, intron-chain fast proposals retain their own junction
  geometry and must agree with their CIGAR reference skips. A fast proposal
  also permits independent generic spliced candidate generation before common
  rescoring; it no longer returns early and hides supported alternatives.
  Scoring constants, coverage fraction, model and evaluator are unchanged.

- In `1.1.1-dev.4`, read-coverage classes protect full candidates from being
  erased by shorter high-scoring explanations. Primary selection continues
  past ineligible score leaders, requires positive score and informative query
  coverage, and keeps different-locus short evidence in the ambiguity margin.
  Deleted reference columns no longer fill read coverage. Mate-window and
  hyperediting arbitration apply the same rules; scoring constants and the
  existing coverage fraction are unchanged. See
  `docs/candidate-eligibility-contract.md`.
- In `1.1.1-dev.3`, split/span rescue proposals are compared with their primary
  using strict read-placement scores. A rewrite must improve sequence evidence
  and exceed the existing split cost; a junction prior alone cannot displace an
  equally good contiguous match. Segmented chains also retain the continuous
  interpretation at each collinear segment even when splice construction
  succeeds. This prevents a repeat-derived splice from hiding the better
  contiguous placement before arbitration. Scoring constants and model weights
  are unchanged; candidate generation and rescue eligibility change.
- In `1.1.1-dev.2`, splice-tail and run-through rescue require both sides to
  meet at the declared query and reference boundaries before accepting a
  stitched CIGAR. Invalid local-extension seams fall back to the original
  placement; valid outer soft clips and indels are retained. Right-tail
  reextension also updates the alignment start and uses that start for the
  short-tail mismatch gate. See `docs/cigar-rescue-contract.md`.
- Score-stage likelihood quadrature no longer stalls on ultra-high-coverage sites. Sites with more than 63 observations were integrated with an adaptive Simpson rule whose refinement cost grows unpredictably on the sharply peaked likelihoods produced by extreme coverage (millions of observations at a single position), leaving affected runs effectively stuck in the score stage. All sites now use fixed 32-point Gauss-Legendre quadrature, so per-site cost scales linearly with observation count and extreme-coverage sites finish in well under a second. Site scores and VCF output are bit-identical to 1.1.0 on full-depth human RNA-seq runs.

## [1.1.0] - 2026-09-29

### Changed

- Candidate discovery and scoring use a unified reference-aware molecular evidence model on the CPU.
- The model bundle contains `model.txt`, `freeze.json`, and `alu.bed`.
- Discovery includes A-to-G and T-to-C observations with MAPQ at least 20 and base quality at least 13; the frozen classifier threshold determines the site filter.
- VCF output reports `ESP_SCORE`, `VAF`, `DEPTH`, `AC`, and `STRAND`, with `PASS` or `LOW_SCORE`. Scores are classifier outputs, not calibrated probabilities.
- Scientific result tables are deterministic across worker counts.
- Resumable runs use metadata schema version 2 and require compatible cached artifacts.

## [1.0.8] - 2026-09-23

### Fixed

- Report generation is now deterministic for multi-copy gene annotation: a gene name present on more than one contig (e.g. PAR genes on chrX/chrY) is attributed to the last GTF row in file order, as documented, instead of whichever copy happened to win a hash-map iteration. Repeated report generation over the same run directory is byte-identical; site tables and scores are unaffected.
- Rescue resume is crash-safe: the rescue completion key is written only after the merged BAM is durable, and `resume` verifies that the sealed rescued count, the rescued sidecar and the `RE:Z:collapsed` records in the BAM all agree before accepting a rescue stage as done. A run interrupted around the merge re-runs the rescue to a consistent state instead of continuing over a half-merged BAM. Resume semantics change: a previously "complete" rescue whose artifacts disagree is now rebuilt (rescued reads are recovered from the unmapped pool or from their own earlier `RE:Z:collapsed` records, with duplicate pool entries dropped), and a stale `RE:Z:collapsed` record whose read no longer verifies is removed at the next merge.

## [1.0.7] - 2026-09-22

### Fixed

- Hyperedited site aggregation now normalizes each rescued read to its genomic polarity before grouping: mate pairs rescued by opposite passes over the same event merge into one cluster and one site row (previously split into AG/TC groups with halved support and half the VCF rows written with the wrong REF/ALT polarity). A tripwire cross-checks every kept site against the reference, reports the outcome in `hyper/polarity_report.json`, and prints a warning listing any inconsistent coordinates (the run still completes); a reference that cannot be read is a hard error.
- Duplicate hyper sites are rejected instead of silently merged or overwritten: the hyper stage refuses to write a sites table with a duplicated (chrom, pos) key, and readers of `hyper/sites.tsv` fail on duplicates, unknown directions, or malformed fields.
- Rescued records written to the merged BAM now carry offset-free Phred quality values; the second `Record::set` overwrote the converted values with raw ASCII bytes, inflating every rescued read's quality by +33.
- `scan` writes `"snr": null` instead of invalid `inf` JSON when a control channel is empty (crashed the downstream verify/score parse).
- Strict input handling: a truncated or malformed `map/unmapped.fq.gz` record and a malformed `map/rescued.tsv` row are hard errors naming the file, the 1-based line and the reason (previously such content was silently dropped); `hyper/sites.tsv` rows with a wrong column count are rejected; a present-but-unreadable `sites.tsv`, `candidates.verified.bed`, `hyper_scores.tsv` or `species.json` is an explicit error instead of being treated as absent.

### Changed

- Removed the `--snp` option from `run` and `scan` (known-SNP exclusion remains an internal channel; the refs-discovery dbSNP default silently altered calling protocols). `resume` now also ignores a known-SNP setting frozen in an older run.json and prints a warning, so no entry point can re-enable the channel; scan-stage artifacts that may have been produced under the filter are rebuilt from scan onward.

## [1.0.6] - 2026-09-14

### Changed

- Updated site scoring and added scores for rescued editing sites.
- Scoring now runs on the CPU; `--device` and `--caduceus` are no longer supported.

## [1.0.5] - 2026-09-11

### Changed

- Improved alignment, candidate calling and hyperediting detection.

## [1.0.4] - 2026-08-31

### Changed

- Improved read rescue, site scoring and report navigation.

## [1.0.3] - 2026-08-29

### Added

- Added hybrid-reference support for human knock-in mouse samples.

## [1.0.2] - 2026-08-27

### Added

- Added resumable runs, per-sample output directories, and update and delete commands.

## [1.0.1] - 2026-08-27

### Added

- Added reference setup, HTML reports and hyperedited-read recovery.

## [1.0.0] - 2026-08-24

### Added

- Initial release with read processing, alignment, editing-site detection and scoring.
