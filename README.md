# ESPERANTO

[![CI](https://github.com/intell-biol-lab/esperanto/actions/workflows/ci.yml/badge.svg)](https://github.com/intell-biol-lab/esperanto/actions/workflows/ci.yml)
[![Development](https://img.shields.io/badge/version-1.1.1--dev.1-orange.svg)](CHANGELOG.md)
![Language: Rust](https://img.shields.io/badge/language-Rust-orange.svg)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

ESPERANTO identifies and scores candidate A-to-I RNA editing sites from RNA sequencing reads using read-level and reference-aware molecular evidence.

Development is maintained at **[intell-biol-lab/esperanto](https://github.com/intell-biol-lab/esperanto)**,
forked from [Zhazhupai603/esperanto](https://github.com/Zhazhupai603/esperanto).
This organization development line has no binary release yet. See the
[development contract](docs/development.md), [contribution guide](CONTRIBUTING.md)
and [baseline audit](docs/provenance/2026-10-29-adoption.json).

**Known qualification failure:** inherited mapper output can contain internal
soft clips in primary CIGARs. The initial human development pilot failed strict
SAM interface qualification. The organization foundation does not repair this
issue or establish benchmark superiority. See [known issues](docs/known-issues.md).

## Installation

Build from source with the pinned Rust toolchain in `rust-toolchain.toml`.
Requires Linux x86-64, a C/C++ toolchain, CMake, pkg-config, libclang, and the
OpenSSL, zlib, bzip2 and xz development libraries. No GPU is required.

```sh
git clone https://github.com/intell-biol-lab/esperanto.git
cd esperanto
# For an analysis, check out its recorded, reviewed commit before building.
cargo build --release --locked -p esperanto-cli
./target/release/esperanto --version
export PATH="$PWD/target/release:$PATH"
# Keep model selection explicit, especially on hosts with an existing install.
export ESPERANTO_BUNDLE="$PWD/bundle/human/esperanto-model-v1.1.0/rust"
```

The inherited Cargo configuration targets x86-64-v3 (AVX2). For a portable
x86-64 executable, set `RUSTFLAGS="-C target-cpu=x86-64"` when building; CI uses
that portable target. Record this choice with each experiment.

`--version` reports the organization repository, exact source commit and build
worktree state; `-V` reports the short software version. Development self-update
is disabled. `install.sh` explains source installation and exits without
installing an upstream binary.

## Quick start

Prepare GRCh38 reference files and indexes once:

```sh
export ESPERANTO_REFS="$PWD/refs"
esperanto setup
export ESP_INDEX="$ESPERANTO_REFS/GRCh38.primary_assembly.genome.paidx"
```

Run paired-end FASTQ input:

```sh
esperanto run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --index "$ESP_INDEX" --sample sample --out results
```

Run a coordinate-sorted and indexed BAM:

```sh
esperanto run --bam sample.sorted.bam --sample sample --out results
```

Use `--sites sites.tsv` with a BAM to score a supplied list. Each TSV row contains a chromosome and a 1-based position, separated by a tab, without a header.

## Commands

| Command | Purpose |
|---|---|
| `setup` | Prepare reference files and indexes |
| `index` | Build alignment indexes from FASTA and optional GTF |
| `run` | Analyze FASTQ or BAM input |
| `resume` | Continue an interrupted run |
| `qc` | Trim and filter FASTQ reads |
| `map` | Align reads to the reference |
| `scan` | Discover candidate editing sites |
| `score` | Score a supplied site list |
| `report` | Regenerate an HTML report |
| `update` | Update the installation |
| `delete` | Remove installed reference or model files |

Run `esperanto --help` or `esperanto run --help` for command options.

## Inputs

- Reference: GRCh38/hg38 FASTA with a matching `.fai`.
- FASTQ: plain or gzip-compressed single-end or paired-end reads; FASTQ runs require the ESPERANTO alignment index.
- BAM: coordinate-sorted and indexed alignments to the selected reference.
- Annotation: an optional matching GTF for splice-junction annotation and HTML reports.
- Model: a complete unified bundle containing `model.txt`, `freeze.json`, and `alu.bed`.

## Outputs

For `--out results --sample sample`, files are written under `results/sample/`.

| File | Contents |
|---|---|
| `sites.vcf` | Candidate sites, classifier scores, filters, and base counts |
| `sample.report.html` | Interactive gene and site report when a GTF is available |
| `qc/*.fq.gz` | Cleaned reads for FASTQ input |
| `map/sorted.bam` | Sorted alignments for FASTQ input |
| `scan/candidates.bed` | Observed A-to-G and T-to-C sites with read counts |
| `score/scores.tsv` | Chromosome, position, and classifier score |
| `run.json` | Input fingerprints and resolved run parameters |

## Behavior and reproducibility

Candidate discovery includes reference A or T positions with at least one editing-consistent observation after MAPQ 20 and base-quality 13 filtering. The unified model supplies the final site decision. Scores are classifier outputs, not calibrated probabilities. `PASS` indicates that `ESP_SCORE` meets the frozen model threshold; `LOW_SCORE` indicates that it does not.

Candidate tables, scores, and VCF records are deterministic for identical inputs and model files, including across worker counts. Run metadata and QC reports contain timing information. `resume` validates saved inputs and stage artifacts before continuing. It accepts run metadata schema version 2; incompatible run directories must be started again in a new output directory.

See [CHANGELOG.md](CHANGELOG.md) for release updates and [technical specifications](docs/specs/crates/flow.md) for file and pipeline contracts.

## Citation

Cite the exact organization commit or organization release used in the analysis,
the inherited upstream software, and the model SHA-256 / threshold separately.
Historical upstream releases remain at
[Zhazhupai603/esperanto](https://github.com/Zhazhupai603/esperanto/releases);
organization release policy is in [docs/development.md](docs/development.md).

## License

[MIT](LICENSE).
