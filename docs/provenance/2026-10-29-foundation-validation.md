# Foundation validation

Validation date: 2026-10-29 (Asia/Shanghai).

## Baseline and change scope

- Inherited software: 80 workspace tests pass, plus one child-process execution
  of the scalar-dispatch test.
- Baseline formatting did not pass rustfmt. Commit
  `741f43a140be296169fe728304236a480660850f` normalizes formatting separately.
- The foundation does not subsequently change mapper, engine, scoring or bundled
  model code/bytes relative to that formatting commit.
- Cargo.lock only changes the ten local workspace package versions;
  third-party package entries and their checksums remain unchanged.
- The [adoption manifest](2026-10-29-adoption.json) records the inherited
  source/model differences and exact hashes.

## Regression tests

The two new `crates/cli/tests/development_identity.rs` integration tests were
run against the old implementation first:

- The old `--version` lacked a development identity and source commit.
- The old updater attempted network access instead of refusing a development
  self-replacement.

Both failed for those expected reasons. Both pass after implementation.

## Local gate

```sh
CARGO_HOME=/data/merlin/.cache/esperanto-cargo CARGO_BUILD_JOBS=8 sh scripts/check.sh
```

The exact path above is the adoption host's writable dependency cache, not a
runtime dependency. Another developer can use their normal Cargo cache.

Results: formatting passes; 82 workspace tests pass, plus one scalar child-test
execution; locked CLI release build succeeds; all three bundled model manifest
entries pass; `--version` and `--help` succeed.

The local gate used Rust 1.94.0, Cargo 1.94.0 and the inherited x86-64-v3 target.
GitHub Actions separately uses the portable x86-64 target.

`sh -n scripts/check.sh install.sh` passes. The installer exits 1 with source
installation guidance. It performs no download or installation.

## Additional build-provenance checks

The actual build.rs was compiled with rustc and executed in temporary real Git
repositories and source directories:

| Case | Observed identity |
| --- | --- |
| Clean checkout | Exact commit, clean |
| Modified tracked source | Same commit, dirty |
| New commit | New exact commit, clean |
| Linked Git worktree | Exact commit, clean |
| Source archive without Git | Unknown commit and worktree |
| Source archive inside an unrelated repository | Unknown commit and worktree |

These checks exercise the build script with actual Git state; they do not
establish bit-identical builds or scientific mapping qualification.

## Remote integration

Foundation PR: https://github.com/intell-biol-lab/esperanto/pull/1

The GitHub API confirmed master requires a PR and up-to-date `verify` status,
including for admins; force pushes and deletion are disabled and conversations
must be resolved. Additional reviewer count is zero for the initial
single-maintainer workflow.

Remote CI is inspected on the current PR head after pushing. Its final result
and exact run/commit identities are retained in the research workspace handoff.
A draft PR is an integration review point, not a stable release.
