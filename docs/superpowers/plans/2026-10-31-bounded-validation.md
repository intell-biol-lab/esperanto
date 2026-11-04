# Bounded validation and cohort evaluation implementation plan

> **For agentic workers:** Use executing-plans to implement the tasks sequentially.

**Goal:** Measure bounded complete-content input hashing and batch amortization,
then independently evaluate candidate-union recount with new multi-sample truth.

**Architecture:** Parallelize independent file digests before any sample worker
starts. Preserve deterministic identities and existing serial output validation.
Keep cohort generation, independent truth audit and evaluation in the research
workspace; they do not change ESP discovery, counting or classifier decisions.

**Tech Stack:** Rust 1.94.0, standard scoped threads, SHA256, Python/pysam.

**Spec:** `docs/specs/crates/batch-joint.md` and the research protocol
`configs/esp-validation-cohort-protocol-20261031.json`.

## Global constraints

- Base: dev.39 `cd6e2d56203e6f21fae501651d0368eacba875cb`.
- Frozen dev.36 and adoption lock remain unchanged. No merge, release or tag.
- Hash complete content; no mtime cache, sampling or skipped files.
- Use at most min(requested validation threads, total threads, unique files).
- Each active hashing worker holds one 1 MiB buffer and one input file.
- Parallel validation finishes before sample workers. Output hashing stays serial.
- No changes to mapping, search, MAPQ, model, discovery or counting thresholds.
- Existing isolated TEST is never opened or evaluated.
- All new evidence goes in `esp-validation-cohort-20261031` roots.

### Task 1: bounded deterministic digests

**Files:** `crates/cli/src/evidence.rs`, `crates/cli/src/batch.rs`,
`crates/cli/tests/batch.rs`.

- [ ] Add tests for serial/parallel identity equality, aliases, content mutation,
  worker cap, invalid cap and filesystem failures. Run RED.
- [ ] Implement scoped worker hashing with deterministic path/result ordering.
- [ ] Add `batch --validation-threads` and per-invocation phase receipts.
  Operational cap changes must not invalidate otherwise identical resume input.
- [ ] Run focused tests, update version/changelog/spec; commit on topic branch.

### Task 2: complete-cost performance evidence

**Files:** new research experiment harness, receipts and report.

- [ ] Build clean committed candidate; retain binary and source provenance.
- [ ] Alternate caps 1/2/4/8 on identical real human resources, three repetitions.
  Compare digest maps independently; record wall, CPU, RSS and I/O.
- [ ] Compare serial versus bounded validation at batch sizes 1/8/32, repeated,
  identical workers, same full verification. Record first-worker delay, total
  wall, output verification and per-sample cost. Verify output equality.
- [ ] Recheck strict mapping and independent audit on the unchanged development
  subsets used by dev.39. Keep performance and accuracy conclusions separate.

### Task 3: independent multi-sample development truth

**Files:** new generator, separate audit and evaluator in research workspace.

- [ ] Freeze generator seed, sample metadata, event classes and decision policies
  before generation. Use new human reference loci and molecular draws.
- [ ] Include shared, condition-specific, private/rare and weak events, depth
  variation, absent coverage, random errors, shared artifact and donor SNP controls.
- [ ] Audit BAM observations and truth from raw reads without importing generator
  or ESP helpers. Reject failed truth before evaluation.
- [ ] Compare unchanged independent discovery against union recount using the
  same predeclared sample-level evidence rule. Report sample-site and site-level
  recall, empirical FDP (an FDR estimate), rare retention and missing states.
- [ ] Evaluate donor removal and sample-order invariance. Never interpret evidence
  states as calibrated editing calls; no claim of validated biological FDR.

### Task 4: review and closure

- [ ] Run `sh scripts/check.sh`, verify old evidence protections.
- [ ] Publish organization draft PR with exact commit and local/remote checks.
- [ ] Record all commands, input/resource/source/binary hashes, limitations and
  result files; externally seal new evidence. Keep dev.36 frozen.
