#!/usr/bin/env sh
set -eu

cat >&2 <<'EOF'
ESPERANTO organization development currently requires a source build.
No organization binary release has been published.

  git clone https://github.com/intell-biol-lab/esperanto.git
  cd esperanto
  # Check out the exact reviewed commit required by your analysis.
  cargo build --release --locked -p esperanto-cli
  ./target/release/esperanto --version

See https://github.com/intell-biol-lab/esperanto for development status.
EOF
exit 1
