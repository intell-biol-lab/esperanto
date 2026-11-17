# Sample batches and independent joint recount

Accepted direction: resource-aware sample-sheet execution first; joint candidate
union and per-sample recount second; cross-sample alignment remains a separate
experimental gate. Parent is frozen dev.36, commit
`28d32aed60429b75664db4a2225a964ffead4f9c`. No mapper, model, threshold, adoption
lock, frozen evidence or held-out TEST changes.

## Executable contract

`esperanto batch --samples samples.tsv --out results --total-threads 32
--threads-per-sample 4 --memory-gib 96 --memory-per-sample-gib 12`
launches independent child processes, with bounded concurrency
`min(total_threads / threads_per_sample, memory_gib / memory_per_sample_gib,
jobs_if_supplied)`. Memory is a reservation estimate, not an OS-enforced limit.
Reject zero resources and impossible reservations. Do not auto-tune from truth.
Default mode `run`; `--mode map` executes the unchanged mapping entry. Shared
reference/model paths must be explicit and are content-fingerprinted once per
batch; child commands include resolved absolute paths. No stage-level overlap
heuristic or resident worker is introduced without measurements.

TSV header `sample`, with `r1,r2,bam,sites,candidates` and optional
`batch,condition,donor,replicate,strandedness`. Unknown/duplicate columns,
duplicate/unsafe sample names, malformed rows, mixed FASTQ/BAM entries, unequal
lane counts, missing files and duplicate physical inputs across samples fail
before launching children. Comma-delimited lanes are for one sample only.
Relative paths resolve against the sample sheet, not the working directory.
Optional metadata is retained; technical batch never implies biological grouping.

One child per sample; no shell interpolation. One isolated output directory and
attempt log per sample. Atomic batch manifest/status records include commands,
sample metadata, parent/child executable hashes, inherited ESP control variables,
all input and resource file SHA-256 identities, configured budgets and timing.
A batch lock excludes concurrent writers; interrupted lock recovery is explicit.
Nonzero child exit is a failed sample, other samples continue, batch exits
nonzero. Resume requires exact pinned identity, validates successful output
hashes, skips intact successes and restarts failed/incomplete samples.
Full-pipeline children with run.json use existing `resume`; map retries require
an empty output path, with previous attempts preserved. Logs retain every attempt.
No successful output may be silently reused after input or configuration drift.

`esperanto joint --samples mapped.tsv --fasta ref.fa --out joint` requires BAMs.
An optional `candidates` column identifies each sample's candidate BED (0-based,
one-base intervals, e.g. existing scan output); otherwise use unchanged discovery.
Build a deterministic candidate union and recount each sample independently
using existing unified molecular observation filters (MAPQ20, BQ13, RG/template
overlap handling). Record candidate provenance and sample metadata separately.
Validate BAM/reference contig lengths; absent BAM contigs are explicit missing
evidence. Reject malformed/out-of-range/non-A-or-T candidates. Never interpret
no coverage as evidence of no editing.

Long-form evidence matrix: chrom, 1-based pos, REF, expected ALT (A>G/T>C),
sample, A/C/G/T counts, informative depth, ALT count/fraction and state.
States: CONTIG_ABSENT; INSUFFICIENT_COVERAGE (< configurable minimum depth,
default 10); ALT_OBSERVED; NO_ALT_OBSERVED. Zero depth has NA fraction, not 0.
Conflicting equal-quality overlapping bases do not contribute informative depth.
Discovery membership is not an editing call; model scores, differential tests,
group pooling and cross-sample MAPQ changes are outside this first implementation.
Output is published only after all recounts succeed, with completion receipt last.

## Acceptance

Behavioral RED tests precede implementation. Test strict sheet parsing and
resource limits; child scheduling, failure isolation, safe argv, output/log
separation, drift rejection and resume; synthetic joint union, overlapping mates,
RG identity, quality/flag filters, A>G/T>C, absent/zero/low-depth evidence,
sample-private candidates, common artifacts, contig mismatch and order stability.
Integration must compare actual candidate batch output with standalone frozen
dev.36, plus independent SAM qualification on named development inputs.
New experiments use a new research directory with exact identities; no claim
of speedup or scientific accuracy gain without a matched experiment.
