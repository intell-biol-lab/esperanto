# Development, provenance and release contract

## Repository ownership

`https://github.com/intell-biol-lab/esperanto` is the continuing development
repository. The upstream repository remains
`https://github.com/Zhazhupai603/esperanto`. Keep `origin` and `upstream` distinct.
Use `master` for integration and topic branches / PRs for changes.

The organization foundation starts from
`3d3c3adfae57f8f8164d1dc5dd5f988f8a18c865`. The earlier frozen publication
source is `d37fda36c16a78d0c403e0b287bcae3a14702192`; all 105 archived source
files match that commit. The fork start differs from those files only in
CHANGELOG. It has additional lockfile and bundled model files.

The inherited `v1.1.0` release points to
`8f955efea1cb949296cf17b5149ac383e1b08a7a`. Its model weights and threshold
differ from the later publication / fork model. Neither the old software version
nor the model folder name is a complete scientific identity.

See the machine-readable [adoption audit](provenance/2026-10-29-adoption.json).

## Build identity

The initial organization version is `1.1.1-dev.1`; it is unreleased. Future
development increments the prerelease identity as needed; every analysis also
pins a full Git commit. Keep historical upstream tags unchanged.

Use Rust 1.94.0 and `--locked`. The workspace package version changes in
Cargo.lock together with Cargo.toml, without upgrading third-party dependencies.
Do not silently regenerate the dependency graph.

`esperanto --version` embeds repository, commit and worktree state at build time.
`-V` remains the short package version. Dirty binaries require a saved patch /
source snapshot and are not release artifacts. When Git is unavailable or the
source is a tarball, commit and worktree state are `unknown`; supply external
archive provenance rather than attributing the build to a surrounding repo.

Build scripts watch source and Git metadata so a changed commit / index updates
the identity on the next build. Changing source after building does not change
an already-built binary. Save the executable SHA-256 with each run.

The local default CPU target is inherited x86-64-v3. CI overrides it with
`RUSTFLAGS="-C target-cpu=x86-64"` and exercises runtime SIMD/scalar dispatch
through the existing tests. A build hash is specific to its toolchain, platform
and flags; this foundation does not claim bit-reproducible binaries.

## Model and reference identity

The inherited bundle is `bundle/human/esperanto-model-v1.1.0/rust`:

- model.txt SHA-256:
  `18899a8b21462528cd388e70cd096f1eba7705c48fefaf3be3c407a4a21245d7`
- classifier threshold: `0.9173483910126397`
- Verify with `sha256sum -c manifest.sha256` from its parent folder.

The publication and fork have identical model.txt bytes and threshold, but their
freeze.json metadata schemas and hashes differ. Do not claim identical complete
bundles. The adoption audit records both freeze.json hashes without importing
the publication's training partitions or local metadata.

Set `ESPERANTO_BUNDLE` to the intended checked-out bundle, or pass the explicit
bundle option to scoring / full-pipeline commands. Record the resolved model
directory and all three model file hashes in experiment provenance.

Inherited `setup` reference defaults are GENCODE v48 with its standard decoy
configuration; the human mapping pilot supplied GENCODE v50 and its index
explicitly. Do not use setup defaults to infer the pilot's exact reference.

If no local bundle resolves, inherited setup still has an upstream `v1.1.0`
tarball fallback. That release's model differs from this fork's model.
Source builds include a bundle and should explicitly select it. Auditing /
replacing that fallback belongs to the release-hardening gate below; it must
not be treated as an equivalent model.

## Integration and CI

`sh scripts/check.sh` runs formatting, workspace tests, a locked CLI release
build, model manifest checks and CLI version/help checks. CI uses the pinned
toolchain, a portable x86-64 target and read-only GitHub permissions.
The required status context is `verify`.

Master protection is configured through GitHub, not by the files alone:
require PRs and the up-to-date `verify` check, apply it to admins, and block force
pushes/deletion. The single-maintainer configuration requires zero extra review
approvals. Add reviewers / CODEOWNERS when actual maintainers are designated.

This CI gate verifies software consistency; the known CIGAR qualification
failure is tracked separately and remains a blocker for a qualified release.

## Releases

There is no organization binary release at adoption. Source installation is
the supported path. install.sh exits with that explanation instead of installing
an upstream binary. Development `update` and `update --check` fail before network
access, including with `--yes`. The stable updater's URLs now name this
organization, but the inherited updater is not certified for an organization
stable release.

Before a first organization stable release:

1. Repair the known CIGAR defect with reproducing tests and pass strict
   independent qualification on both human development seeds.
2. Review mapping deltas, read / molecule denominators, and model / scoring
   compatibility. Do not claim general superiority from the small pilot.
3. Audit installer and updater checksums, archive extraction, atomic replacement,
   semantic version comparison (including prereleases) and model compatibility.
   Replace or explicitly remove the differing upstream model fallback.
4. Verify clean-checkout and source-tarball builds and run identities.
5. Document model, reference, CLI and resume compatibility. The organization
   version bump intentionally prevents silent resume of old versioned runs.
6. Require CI, a release note, immutable tag/commit, portable build description,
   binary / archive / model checksums and a model card. Obtain release
   authorization before tagging or publishing.

Existing upstream Git tags are historical references; a forked tag is not a
new organization software release.
