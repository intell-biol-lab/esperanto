# Organization Foundation Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to execute this plan
> sequentially. Do not dispatch subagents without user authorization.

**Goal:** Make the organization fork a reproducible, traceable development home.

**Architecture:** Preserve upstream history, isolate the initial infrastructure
change on a topic branch, and establish one CI gate before protecting master.
Keep mapping semantics and the inherited model fixed.

**Tech Stack:** Rust 1.94.0, Cargo, GitHub Actions, GitHub pull requests.

**Spec:** `docs/superpowers/specs/2026-10-29-organization-foundation.md`

## Global constraints

- Source baseline: `3d3c3adfae57f8f8164d1dc5dd5f988f8a18c865`.
- Rust: `1.94.0`; organization development version: `1.1.1-dev.1`.
- Preserve Cargo.lock's dependency graph, inherited release tags and model files.
- No mapper repair, reference change, retraining or stable release.

## Tasks

- [x] Audit source, release, lockfile and model identities; run inherited tests.
- [x] Normalize inherited formatting in its own commit; establish formatting gate.
- [x] Add failing CLI integration tests for development version identity,
  commit reporting and offline refusal of self-update.
- [x] Add development version, build provenance and self-update guard.
- [x] Add pinned toolchain, CI workflow, contribution / release contracts,
  pull request and issue templates, baseline identity manifest and fork README.
- [x] Run workspace tests, formatting, locked release build and CLI smoke checks;
  verify model files and mapping behavior source were not semantically changed.
- [x] Commit and push the topic branch, open a draft PR and inspect remote CI.
- [x] Configure master protection using the observed CI check name and verify
  settings. Record the verified result in the research workspace handoff.
