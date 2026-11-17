# ESPERANTO maintenance instructions

The canonical development repository is `intell-biol-lab/esperanto`.
Preserve upstream attribution, Git history, LICENSE and historical tags.

- Work on a topic branch and use a pull request into `master`.
- Read `CONTRIBUTING.md`, `docs/development.md` and `docs/known-issues.md`.
- Use pinned Rust and tracked Cargo.lock. `sh scripts/check.sh` is the CI gate.
- Keep mapping, scoring, model and evaluation changes separately reviewable.
- For fixes, add a reproducing regression test before changing behavior.
- Never make invalid CIGAR output pass by weakening the independent scorer,
  silently rewriting internal soft clips, or excluding failed records from
  scientific denominators.
- Do not mutate frozen publication artifacts or existing accepted truth/run
  directories. Put new experiments in new directories and record exact inputs.
- Preserve model weights and thresholds unless a model change is explicitly
  part of the task. Record model identity separately from software version.
- Before claiming qualification or performance, run the independent evaluator
  on the named inputs. Unit tests and a successful build are insufficient.
- Capture commit, dirty state, toolchain, build flags, binary hash, model hashes,
  reference/index identity, commands and result files for scientific runs.
- Do not commit private human reads, credentials, generated indexes or outputs.
- Do not create a release/tag or merge a PR unless authorized by the task.
