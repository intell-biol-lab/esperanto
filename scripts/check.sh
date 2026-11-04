#!/usr/bin/env sh
# The same gate used locally and by GitHub Actions. Run from any directory.
set -eu
cd "$(dirname "$0")/.."

cargo fmt --all -- --check
cargo test --workspace --locked
cargo build --release --locked -p esperanto-cli
(
    cd bundle/human/esperanto-model-v1.1.0
    sha256sum -c manifest.sha256
)
./target/release/esperanto --version
./target/release/esperanto --help > /dev/null
sha256sum target/release/esperanto
