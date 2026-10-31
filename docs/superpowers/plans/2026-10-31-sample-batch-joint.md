# Sample batch and joint evidence implementation plan

> **For agentic workers:** use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** deliver bounded sample-sheet execution and an independent per-sample
candidate recount interface on frozen dev.36.

**Architecture:** shared strict TSV parser in flow; CLI batch coordinator launches
unchanged child entry points; score owns candidate union/recount; CLI joint only
wires paths/provenance. Existing mapper/scorer model behavior stays unchanged.

**Tech Stack:** pinned Rust 1.94, existing serde/serde_json/sha2/rust-htslib.

**Spec:** `docs/superpowers/specs/2026-10-31-sample-batch-joint.md`.

## Global constraints

No mapper/model/evaluator changes, no adoption-lock update, no held-out TEST.
Preserve all sealed evidence. Work on topic branch and organization draft PR.

### Task 1: sheet and execution contract

- [ ] Add `crates/flow/src/samples.rs` parser tests: duplicate names, path traversal,
  relative lanes, mixed entries, malformed metadata and aliased inputs.
- [ ] Run new tests against an unimplemented parser and retain RED output.
- [ ] Implement `read_samples(path) -> Result<Vec<Sample>>`, preserving metadata.
- [ ] Add `crates/cli/src/batch.rs`, resource planning and mock-child integration
  tests for concurrency, argv isolation, failure and content-verified resume.
- [ ] Implement command, identity manifest, lock, attempt logs and status writes.
- [ ] Verify parser and batch tests, then commit this independently reviewable task.

### Task 2: independent recount

- [ ] Add scientific tests in `crates/score/src/joint.rs` using synthetic indexed
  BAM/FASTA with overlapping mates and sample-specific candidates.
- [ ] Demonstrate RED for candidate union and three coverage states.
- [ ] Add read-only `unified::count_sites` using the existing observation reader.
- [ ] Implement `joint::recount` and CLI `joint`, deterministic outputs, metadata
  and completion receipt; verify validation and sample-add/remove invariance.
- [ ] Verify score/joint tests and commit separately.

### Task 3: development qualification and review

- [ ] Run complete `sh scripts/check.sh`, with new CLI integration tests.
- [ ] Save exact clean build/source/model/input identities in a fresh experiment.
- [ ] Compare frozen dev.36 standalone with candidate batch on named dev reads,
  including all SAM fields/tags and independent mapping qualification.
- [ ] Run joint against independent synthetic expected counts and audit BAM hashes.
- [ ] Write result report, usage docs and remaining experimental gate.
- [ ] Push topic branch and create organization draft PR; record CI status.

Cross-sample alignment and inferential joint statistics require new qualified
multi-sample truth before implementation/adoption; this task supplies their
input/evidence contract and cannot claim that scientific gate has passed.
