# Batch execution and joint evidence

Status: development, on top of frozen dev.36. Neither command changes mapping
candidate selection, scores, MAPQ, model weights or per-sample calling thresholds.

## Sample sheets

UTF-8 TSV, header required, one row per sample. Paths resolve relative to the sheet.
Comments starting `#` and empty lines are ignored. Supported columns:

| Column | Meaning |
| --- | --- |
| sample | Unique identifier, starts with alphanumeric; then alphanumeric, `.`, `_`, `-`; at most 128 bytes |
| r1, r2 | Comma-separated lanes of one sample; R2 lane count equals R1 |
| bam | One coordinate-sorted indexed BAM, mutually exclusive with FASTQ |
| sites | `run` only: existing 1-based `chrom<TAB>pos` TSV |
| candidates | `joint` only: existing one-base BED, 0-based start/exclusive end |
| batch, condition, donor, replicate | Separate metadata, retained verbatim; no automatic biological pooling |
| strandedness | empty, unstranded, forward, reverse; metadata only, does not reconfigure the mapper |

Unknown/duplicate headers, bad rows, duplicate names, missing files and shared
FASTQ/BAM physical files (including Unix hard links) fail before workers start.
Paths cannot contain commas/newlines/tabs. Sample names cannot contain spaces.
Candidate files can be shared intentionally; BAM/FASTQ files cannot.

## Batch

```sh
esperanto batch --samples samples.tsv --out results \
  --index human.paidx --fasta human.fa --gtf genes.gtf --l1-bundle human.bndl \
  --bundle bundle/human/esperanto-model-v1.1.0/rust \
  --total-threads 32 --threads-per-sample 4 \
  --memory-gib 96 --memory-per-sample-gib 12
```

`--mode run` (default) launches the existing full pipeline, with isolated
`results/samples/<sample>/`. `--mode map` launches unchanged `map`, with exactly
one R1 and optional R2 per sample. `--stop-after` is run-only.
`--executable` pins an alternative ESP worker, including a frozen baseline;
the coordinator and worker binary identities are both retained.

Concurrency is the minimum of CPU-thread and estimated memory reservations,
optionally capped with `--jobs`. These are scheduling budgets, not CPU affinity
or enforced memory caps. Auxiliary threads and actual memory use still need
deployment measurements; leave memory headroom. No automatic NUMA placement,
stage scheduling, common discovered-junction sharing or resident worker is used.
Immutable mmap pages can already be shared by the operating system.

`--dry-run` hashes inputs/resources and publishes the reviewable plan, without
launching children. Continue with the same command plus `--resume`.
`batch.json` records resolved commands, metadata, all file hashes, selected
ESP environment controls, resource reservations and per-attempt status/timing.
`logs/` retains argv/stdout/stderr for every attempt. Successful child execution
does not substitute for independent mapping qualification.

Failures do not cancel unrelated samples; the coordinator exits nonzero when
any sample failed. `--resume` verifies exact batch identity and all completed
outputs, skips intact successes and retries failed/pending work. Existing full
runs use native `resume` (except stop-after jobs, which restart in a fresh
directory); map/other retries archive prior directories under `failed-attempts/`.
The parent `.lock` excludes concurrent coordinators. After a killed coordinator,
confirm its workers have stopped before explicitly removing a stale lock;
native per-sample run locks may also need that check. Normal coordinator errors
kill and reap its immediate active children.

All content hashing is part of end-to-end batch overhead. It is not hidden from
user-facing timing; it can dominate very small batches with large references.
This version does not introduce an unverified identity cache to avoid that cost.

## Joint recount

```text
sample	bam	candidates	batch	condition	donor
A	../results/samples/A/map/sorted.bam	A.bed	B1	case	P1
B	../results/samples/B/map/sorted.bam	B.bed	B1	control	P2
```

```sh
esperanto joint --samples mapped.tsv --fasta human.fa --out joint \
  --min-depth 10 --threads 4
```

Omit `candidates` to run unchanged per-sample discovery. `sites` is rejected here
to prevent silently confusing 1-based positions with BED coordinates.
BAM contig names/lengths must match the supplied reference; absent BAM contigs
are explicit missing evidence. This validates dictionaries, not reference
sequence identity within old BAM headers; users must supply the correct FASTA.
Use matching sequence hashes/provenance where available.

Outputs:

- `candidates.bed`: deterministic union, including sample-private candidates.
- `candidate_origins.tsv`: discovery membership, not an editing call.
- `evidence.tsv`: long-form per-sample A/C/G/T counts and coverage state.
- `samples.json`, `provenance.json`: metadata, identities and counting policy.
- `completion.json`: published last, includes hashes of all completed outputs.

Counting reuses the existing unified molecular reader: MAPQ >=20, BQ >=13,
exclude unmapped/secondary/supplementary, merge RG/template overlapping bases,
highest-quality observation wins and equal-quality base conflicts are omitted.
The frozen reader retains duplicate-flagged records; this is not UMI deduplication.
Counts are never pooled across samples. A>G and reference-strand T>C are counted;
transcript-strand reassignment is not inferred from sheet metadata.

| State | Meaning |
| --- | --- |
| CONTIG_ABSENT | Candidate contig absent in this BAM; counts/fraction NA |
| INSUFFICIENT_COVERAGE | Informative depth below min-depth; observed counts retained |
| ALT_OBSERVED | Sufficient depth, expected ALT count >0 |
| NO_ALT_OBSERVED | Sufficient depth, expected ALT count =0 |

Zero-depth ALT fraction is NA. These are evidence states, not probability-based
editing calls. No classifier rescore, differential test, group shrinkage or
cross-sample alignment assistance is introduced. Adding another sample can
expand the union, but cannot change existing per-sample counts at an old site.
Changing discovery membership/group metadata never supplies ALT observations.

## Next scientific gate

Build independent development and sealed test multi-sample truth with shared,
private, condition-specific, low-depth and common-artifact events; include
different donors and technical batches. Compare independent discovery with union
recount, quantify recall/FDR per class and test sample-order/add/remove stability.
Fit/test splits must isolate related molecules/donors; no TEST-derived candidate
library or parameters may enter development. Only then qualify replicate-aware
statistics and a separate frozen cross-sample splice proposal experiment.
