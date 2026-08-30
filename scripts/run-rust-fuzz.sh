#!/usr/bin/env bash
set -euo pipefail

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"
readonly fuzz_root="${repository_root}/tests/fuzz"
readonly fuzz_toolchain='nightly-2026-07-31'
readonly expected_rustc='rustc 1.99.0-nightly (8ab9fdff5 2026-07-30)'
readonly expected_cargo_fuzz='cargo-fuzz 0.13.2'
readonly fuzz_seed=42574257
readonly max_len=4096
readonly targets=(fasta fastq paired_fastq)
readonly requested_output="${1:-}"
readonly seconds_per_target="${BSBIT_FUZZ_SECONDS_PER_TARGET:-30}"

cd -- "${repository_root}"

if [[ $# -gt 1 ]]; then
  printf 'usage: %s [ABSENT_OUTPUT_DIRECTORY]\n' "$0" >&2
  exit 2
fi
if [[ ! ${seconds_per_target} =~ ^[0-9]+$ ]] \
  || (( seconds_per_target < 1 || seconds_per_target > 3600 )); then
  printf 'error: BSBIT_FUZZ_SECONDS_PER_TARGET must be an integer from 1 through 3600\n' >&2
  exit 2
fi
if [[ $(cargo fuzz --version) != "${expected_cargo_fuzz}" ]]; then
  printf 'error: cargo-fuzz must be exactly %s\n' "${expected_cargo_fuzz}" >&2
  exit 1
fi
if [[ $(rustc "+${fuzz_toolchain}" -V) != "${expected_rustc}" ]]; then
  printf 'error: fuzz nightly does not match %s\n' "${expected_rustc}" >&2
  exit 1
fi
for command in cargo c++ git python3 rg rustc sha256sum; do
  command -v "${command}" >/dev/null || {
    printf 'error: missing fuzz tool: %s\n' "${command}" >&2
    exit 1
  }
done

formal_capture=0
if [[ -n ${requested_output} ]]; then
  formal_capture=1
  [[ ! -e ${requested_output} && ! -L ${requested_output} ]] || {
    printf 'error: output already exists: %s\n' "${requested_output}" >&2
    exit 1
  }
  output_parent="$(dirname -- "${requested_output}")"
  [[ -d ${output_parent} ]] || {
    printf 'error: output parent does not exist: %s\n' "${output_parent}" >&2
    exit 1
  }
  git diff --quiet --
  git diff --cached --quiet --
  [[ -z $(git ls-files --others --exclude-standard) ]] || {
    printf 'error: untracked source files exist; capture only from committed source\n' >&2
    exit 1
  }
  for tracked_input in \
    tests/fuzz/Cargo.toml \
    tests/fuzz/Cargo.lock \
    scripts/run-rust-fuzz.sh \
    scripts/summarize-rust-fuzz.py
  do
    git ls-files --error-unmatch -- "${tracked_input}" >/dev/null
  done
  stage_output="$(mktemp -d "${output_parent}/.$(basename -- "${requested_output}").partial.XXXXXX")"
else
  stage_output="$(mktemp -d --tmpdir bsbit-rust-fuzz.XXXXXX)"
fi
readonly formal_capture
readonly stage_output

cleanup() {
  status=$?
  if (( status != 0 && formal_capture == 1 )); then
    printf 'fuzz failure evidence retained at %s\n' "${stage_output}" >&2
  elif [[ -d ${stage_output} ]]; then
    rm -rf -- "${stage_output}"
  fi
  trap - EXIT
  exit "${status}"
}
trap cleanup EXIT

source_manifest() {
  git ls-files -z -- \
    Cargo.toml \
    Cargo.lock \
    crates/bsbit-align \
    crates/bsbit-core \
    crates/bsbit-hts \
    crates/bsbit-index \
    crates/bsbit-io \
    tests/fuzz \
    scripts/run-rust-fuzz.sh \
    scripts/summarize-rust-fuzz.py \
    | sort -z \
    | xargs -0 sha256sum
}

mkdir -p \
  "${stage_output}/artifacts/fasta" \
  "${stage_output}/artifacts/fastq" \
  "${stage_output}/artifacts/paired_fastq" \
  "${stage_output}/corpus/fasta" \
  "${stage_output}/corpus/fastq" \
  "${stage_output}/corpus/paired_fastq" \
  "${stage_output}/logs"

source_manifest > "${stage_output}/source-before.sha256"
git rev-parse HEAD > "${stage_output}/commit.txt"
{
  printf 'cargo_fuzz=%s\n' "$(cargo fuzz --version)"
  printf 'rustc=%s\n' "$(rustc "+${fuzz_toolchain}" -Vv | tr '\n' ';')"
  printf 'cargo=%s\n' "$(cargo "+${fuzz_toolchain}" -V)"
  printf 'cxx=%s\n' "$(c++ --version | sed -n '1p')"
  printf 'kernel=%s\n' "$(uname -a)"
  printf 'seconds_per_target=%s\n' "${seconds_per_target}"
  printf 'seed=%s\n' "${fuzz_seed}"
  printf 'max_len=%s\n' "${max_len}"
  printf 'rustflags=-A dead-code\n'
  printf 'asan_options=detect_leaks=0:halt_on_error=1:abort_on_error=1\n'
} > "${stage_output}/environment.txt"

printf '\001>chr1\nACGTN\n' > "${stage_output}/corpus/fasta/valid"
printf '\037>\n' > "${stage_output}/corpus/fasta/malformed"
printf '\001@read/1\nACGTN\n+\nABCDE\n' > "${stage_output}/corpus/fastq/valid"
printf '\003@truncated\nACGT\n' > "${stage_output}/corpus/fastq/truncated"
python3 -c 'from pathlib import Path; import sys; p=Path(sys.argv[1]); a=b"@pair/1\nACGT\n+\nABCD\n"; b=b"@pair/2\nTGCA\n+\nDCBA\n"; p.write_bytes(bytes((1,7))+len(a).to_bytes(2,"little")+a+b)' \
  "${stage_output}/corpus/paired_fastq/valid"
python3 -c 'from pathlib import Path; import sys; p=Path(sys.argv[1]); a=b"@left/1\nA\n+\nA\n"; b=b"@right/2\nT\n+\nB\n"; p.write_bytes(bytes((2,5))+len(a).to_bytes(2,"little")+a+b)' \
  "${stage_output}/corpus/paired_fastq/name-mismatch"

export CARGO_NET_OFFLINE=true
export ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:abort_on_error=1'
# Rust 1.94's diagnostic renderer crashes under the accepted host's bogus
# 131072x1 WSL PTY while rendering existing dead-code warnings. Fuzz targets
# do not use those diagnostics as an oracle, so suppress only that lint.
export RUSTFLAGS='-A dead-code'
python3 scripts/summarize-rust-fuzz.py self-test
cargo fmt --manifest-path "${fuzz_root}/Cargo.toml" -- --check
cargo clippy --manifest-path "${fuzz_root}/Cargo.toml" --locked --offline --all-targets -- -D warnings

run_fuzzer() {
  local command=$1
  shift
  cargo "+${fuzz_toolchain}" fuzz "${command}" \
    --fuzz-dir "${fuzz_root}" "$@"
}

run_fuzzer build

for target_name in "${targets[@]}"; do
  log_path="${stage_output}/logs/${target_name}.log"
  run_fuzzer run \
    "${target_name}" \
    "${stage_output}/corpus/${target_name}" \
    -- \
    "-max_total_time=${seconds_per_target}" \
    "-max_len=${max_len}" \
    -timeout=5 \
    -rss_limit_mb=2048 \
    "-seed=${fuzz_seed}" \
    -print_final_stats=1 \
    "-artifact_prefix=${stage_output}/artifacts/${target_name}/" \
    > "${log_path}" 2>&1
  if [[ -n $(find "${stage_output}/artifacts/${target_name}" -type f -print -quit) ]]; then
    printf 'error: fuzz target produced a failure artifact: %s\n' "${target_name}" >&2
    exit 1
  fi
done

python3 scripts/summarize-rust-fuzz.py summary \
  --requested-seconds "${seconds_per_target}" \
  --seed "${fuzz_seed}" \
  --max-len "${max_len}" \
  --output "${stage_output}/summary.tsv" \
  "${stage_output}/logs/fasta.log" \
  "${stage_output}/logs/fastq.log" \
  "${stage_output}/logs/paired_fastq.log"

source_manifest > "${stage_output}/source-after.sha256"
cmp -- "${stage_output}/source-before.sha256" "${stage_output}/source-after.sha256"
scripts/check-native-sources.sh > "${stage_output}/native-sources.txt"
printf 'PASS\n' > "${stage_output}/result.txt"
(
  cd -- "${stage_output}"
  find . -type f ! -name SHA256SUMS -printf '%P\0' \
    | sort -z \
    | xargs -0 sha256sum \
    > SHA256SUMS
)

if (( formal_capture == 1 )); then
  mv -- "${stage_output}" "${requested_output}"
  printf 'Rust coverage-guided fuzz artifact published at %s\n' "${requested_output}"
  sed -n '1,5p' "${requested_output}/summary.tsv"
else
  sed -n '1,5p' "${stage_output}/summary.tsv"
fi
