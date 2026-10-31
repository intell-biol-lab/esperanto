# Organization development foundation

Date: 2026-10-29 (Asia/Shanghai).

The maintainer has designated `intell-biol-lab/esperanto` as the continuing
development repository. This work establishes that foundation; the mapper
internal-soft-clip repair remains a separate change.

## Repository and baseline

- Keep `master` as the integration branch and use topic branches and pull requests.
- Use `origin` for the organization repository and `upstream` for
  `Zhazhupai603/esperanto`. Never rewrite inherited history or release tags.
- Start from upstream commit `3d3c3adfae57f8f8164d1dc5dd5f988f8a18c865`.
- Record the publication source commit
  `d37fda36c16a78d0c403e0b287bcae3a14702192` separately from upstream release
  `v1.1.0` (`8f955efea1cb949296cf17b5149ac383e1b08a7a`).
- All 105 files in the frozen source archive match the publication source
  commit. Only CHANGELOG differs between those files and the fork starting
  commit. The published tag and later development tree have different model
  weights and thresholds despite retaining the same model directory name.

## Reproducibility and installation

- Keep Cargo.lock tracked and build/test with `--locked`.
- Pin Rust 1.94.0, including rustfmt and clippy.
- Identify organization development as `1.1.1-dev.1`, with a separate inherited
  model identity. This is not a stable software release.
- Report repository, commit and dirty state in the CLI long version. Source
  archives without Git must explicitly report unknown provenance.
- Do not let development binaries replace themselves with an upstream binary.
  Development self-update fails before any network request; any eventual
  stable updater uses the organization release endpoint.
- Until organization binaries exist, document source installation and make
  the inherited installer fail with an actionable explanation.
- Keep the bundled model and reference defaults unchanged. Document that setup
  still has a pinned upstream model fallback and GENCODE v48 defaults; the human
  mapping pilot used explicitly supplied GENCODE v50 resources.

## Quality and maintenance

- CI runs on pushes and pull requests with read-only repository permissions,
  a portable x86-64 CPU target, the pinned toolchain and locked dependencies.
- Required gate: formatting, workspace tests, CLI release build and provenance /
  development-update smoke checks.
- Normalize inherited formatting in a separate commit so subsequent format
  checks are usable. Do not mix mapper fixes into this work.
- Document contribution, version, model, experiment and release contracts,
  including the known mapping qualification failure.
- Protect master against force pushes/deletion and require the CI check and
  a pull request. Do not require a second person's approval for this currently
  single-maintainer workflow; admins are subject to the gate.
- Keep source data, private human data, model retraining and large indexes out
  of repository changes. A future mapper fix must rerun independent interface
  qualification and retain baseline comparisons.

## Acceptance

The foundation is reviewable as a pushed topic branch and PR, with passing local
checks, a verified remote CI result, an exact source/model baseline manifest and
documented branch settings. No stable release, CIGAR repair or performance claim
is implied by this foundation.
