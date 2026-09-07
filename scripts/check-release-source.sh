#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"

for command in cargo git python3; do
  command -v "${command}" >/dev/null 2>&1 || die "missing command: ${command}"
done

cd -- "${repository_root}"

git diff --check
cargo fmt --all -- --check
cargo fmt --manifest-path tests/fuzz/Cargo.toml --all -- --check
cargo check --locked --manifest-path tests/fuzz/Cargo.toml --all-targets
cargo clippy --locked --workspace --all-targets --all-features --no-deps -- -D warnings
cargo test --locked --workspace --all-targets --all-features
python3 -m unittest discover --start-directory tests/tools --pattern 'test_*.py'

readonly duplicate_dependencies="$(cargo tree --locked --workspace --all-features --duplicates)"
if [[ -n ${duplicate_dependencies} ]]; then
  printf '%s\n' 'error: duplicate dependency versions detected:' >&2
  printf '%s\n' "${duplicate_dependencies}" >&2
  exit 1
fi

printf '%s\n' 'release source checks passed'
