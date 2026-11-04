# Contributing

The organization fork is the continuing development home. `master` is the
integration branch; it is not a claim that all scientific qualification gates
have passed.

## Start work

```sh
git clone https://github.com/intell-biol-lab/esperanto.git
cd esperanto
git remote add upstream https://github.com/Zhazhupai603/esperanto.git
git remote set-url --push upstream DISABLED
git config remote.pushDefault origin
git fetch upstream --tags
git switch -c fix/descriptive-name origin/master
```

Keep upstream history and tags intact. Review upstream commits before importing
them; do not overwrite the organization's integration branch with upstream.
Use `feat/`, `fix/`, `chore/`, `docs/` or `test/` topic branches.

## Changes and checks

Use Rust 1.94.0 and committed Cargo.lock. Install the native dependencies listed
in README. Add regression tests for behavioral changes; run them against the
old implementation to demonstrate the failure.

```sh
cargo fmt --all
sh scripts/check.sh
```

CI builds for portable x86-64; the inherited local default is x86-64-v3.
To reproduce CI, export `RUSTFLAGS="-C target-cpu=x86-64"` first. `ESP_NO_SIMD=1`
can additionally exercise the aligner's scalar dispatch; it does not lower the
CPU requirements of a binary compiled for x86-64-v3.

Keep changes small. Separate formatting, mapper fixes, scoring changes and
model changes. Update CHANGELOG and relevant executable specifications.

## Pull requests

Push to the organization repository and open a PR into `master`. Describe the
trigger, resulting behavior, exact validation and remaining limitations. Link
baseline and output evidence for scientific changes. The `verify` CI job must
pass. Maintainers review and merge; independent human approval is encouraged,
but the initial single-maintainer gate does not require a second reviewer.

Use a merge commit for the foundation PR so its formatting-only commit remains
independently inspectable. Subsequent PRs may use squash merging.
Delete merged topic branches when they no longer support active work.

## Scientific evidence

For each mapping experiment record source commit and dirty state, executable
SHA-256, toolchain / build flags, separate model hashes and threshold,
reference/index hashes, input/truth hashes, full command line and exit status,
and independent evaluator version and outputs.

Keep failed qualification results. A revised mapper gets a new run directory.
Do not overwrite the archived ESP baseline or first pilot outputs.

See [development and release policy](docs/development.md).
