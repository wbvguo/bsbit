#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export CARGO_INCREMENTAL=0

readonly reader_cases=8192
readonly record_cases=8192
readonly header_cases=4096
readonly text_cases=100000

usage() {
  printf 'Usage: %s [OUTPUT_DIRECTORY]\n' "${0##*/}"
  printf 'Default: artifacts/bsbit-release-soak-<UTC timestamp>\n'
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

if (( $# > 1 )); then
  usage >&2
  exit 2
fi
if [[ ${1-} == '--help' || ${1-} == '-h' ]]; then
  usage
  exit 0
fi

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(git -C "${script_dir}/.." rev-parse --show-toplevel)"
cd -- "${repository_root}"

for command in autoreconf cargo cmp cp find gcc git head make mv python3 rg rm \
  rustc sha256sum sort uname xargs; do
  command -v "${command}" >/dev/null 2>&1 \
    || die "required command is unavailable: ${command}"
done
[[ -x /usr/bin/time ]] || die 'required executable is unavailable: /usr/bin/time'

requested=${1:-"artifacts/bsbit-release-soak-$(date -u +%Y%m%dT%H%M%SZ)"}
if [[ ${requested} == /* ]]; then
  final_output=${requested}
else
  final_output="${repository_root}/${requested}"
fi
output_parent=$(dirname -- "${final_output}")
output_base=$(basename -- "${final_output}")
[[ -n ${output_base} && ${output_base} != '.' && ${output_base} != '..' ]] \
  || die "invalid output directory: ${requested}"
mkdir -p -- "${output_parent}"
output_parent="$(cd -- "${output_parent}" && pwd -P)"
final_output="${output_parent}/${output_base}"
[[ ! -e ${final_output} && ! -L ${final_output} ]] \
  || die "refusing to overwrite existing output: ${final_output}"

stage_output="$(mktemp -d -- "${output_parent}/.${output_base}.partial.XXXXXX")"
native_root="$(mktemp -d --tmpdir bsbit-release-soak-native.XXXXXX)"
readonly stage_output final_output native_root
completed=0
cleanup() {
  local status=$?
  if (( completed == 0 )) && [[ -d ${stage_output} ]]; then
    printf '%s\n' "${status}" > "${stage_output}/exit-status.txt"
    if [[ -n ${sanitized_source-} && -f ${sanitized_source}/config.log ]]; then
      cp -- "${sanitized_source}/config.log" \
        "${stage_output}/sanitized-htslib-config.log" || true
    fi
    bash "${script_dir}/check-native-sources.sh" \
      > "${stage_output}/native-sources-after-failure.txt" 2>&1 || true
    printf 'release soak failed; partial evidence retained at %s\n' \
      "${stage_output}" >&2
  fi
  rm -rf -- "${native_root}"
  exit "${status}"
}
trap cleanup EXIT

git diff --quiet --exit-code \
  || die 'tracked worktree changes exist; capture only from committed source'
git diff --cached --quiet --exit-code \
  || die 'staged changes exist; capture only from committed source'
[[ -z $(git ls-files --others --exclude-standard) ]] \
  || die 'untracked source files exist; capture only from committed source'
bash "${script_dir}/check-native-sources.sh" \
  > "${stage_output}/native-sources-before.txt"
readonly -a source_files=(
  Cargo.toml
  Cargo.lock
  crates/bsbit-cli/tests/cli.rs
  crates/bsbit-hts/tests/bam.rs
  crates/bsbit-hts/tests/text_record_fuzz_smoke.rs
  crates/bsbit-hts/htslib-shim/bsbit_hts.h
  crates/bsbit-hts/htslib-shim/bsbit_hts.c
  crates/bsbit-hts/htslib-shim/tests/mutation_smoke.c
  scripts/run-release-soak.sh
)
sha256sum "${source_files[@]}" > "${stage_output}/source-inputs.pre-run.sha256"

sanitized_source="${native_root}/htslib-1.24-sanitized"
mkdir -- "${sanitized_source}"
cp -a -- external/htslib/. "${sanitized_source}/"
rm -rf -- "${sanitized_source}/.git" "${sanitized_source}/htscodecs/.git"
printf '#define HTSCODECS_VERSION_TEXT "1.6.7"\n' \
  > "${sanitized_source}/htscodecs/htscodecs/version.h"
(
  cd -- "${sanitized_source}"
  export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
  export UBSAN_OPTIONS=halt_on_error=1
  autoreconf -i
  CFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \
  LDFLAGS='-fsanitize=address,undefined' \
    ./configure --disable-plugins --disable-libcurl --disable-gcs --disable-s3
  make -j2 lib-static
) > "${stage_output}/sanitized-htslib-build.log" 2>&1
rg -q -- '-fsanitize=address,undefined' "${sanitized_source}/config.mk"

readonly mutation_binary="${native_root}/bsbit_htslib_shim_mutation_soak"
gcc -std=c11 -D_XOPEN_SOURCE=700 \
  -DREADER_CASES="${reader_cases}u" \
  -DRECORD_CASES="${record_cases}u" \
  -DHEADER_CASES="${header_cases}u" \
  -Wall -Wextra -Wpedantic -Wconversion -Werror \
  -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -fno-sanitize-recover=undefined \
  -Icrates/bsbit-hts/htslib-shim -isystem "${sanitized_source}" \
  crates/bsbit-hts/htslib-shim/bsbit_hts.c \
  crates/bsbit-hts/htslib-shim/tests/mutation_smoke.c \
  "${sanitized_source}/libhts.a" -ldeflate -llzma -lbz2 -lz -lm -lpthread \
  -o "${mutation_binary}" \
  > "${stage_output}/native-mutation-build.log" 2>&1
sha256sum "${mutation_binary}" > "${stage_output}/native-mutation-binary.sha256"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
  /usr/bin/time -v -o "${stage_output}/native-mutation-time.txt" \
  "${mutation_binary}" "${native_root}/mutation-artifacts" \
  > "${stage_output}/native-mutation-result.txt"
rg -q "reader_cases=${reader_cases} record_cases=${record_cases} header_cases=${header_cases}" \
  "${stage_output}/native-mutation-result.txt"

BSBIT_TEXT_FUZZ_CASES="${text_cases}" \
  /usr/bin/time -v -o "${stage_output}/rust-text-fuzz-time.txt" \
  cargo test --release --locked -p bsbit-hts --test text_record_fuzz_smoke \
  > "${stage_output}/rust-text-fuzz-result.txt" 2>&1
readonly rust_bam_mutation_test='serialized_bam_mutations_are_bounded_deterministic_and_structurally_checked'
if ! cargo test --release --locked -p bsbit-hts --test bam -- --list \
  | rg -Fqx "${rust_bam_mutation_test}: test"; then
  printf 'error: release-soak BAM mutation test is absent: %s\n' \
    "${rust_bam_mutation_test}" >&2
  exit 1
fi
/usr/bin/time -v -o "${stage_output}/rust-bam-mutation-time.txt" \
  cargo test --release --locked -p bsbit-hts --test bam \
  "${rust_bam_mutation_test}" \
  -- --exact > "${stage_output}/rust-bam-mutation-result.txt" 2>&1
readonly rust_cli_permission_test='input_and_output_permission_failures_publish_nothing'
if ! cargo test --release --locked -p bsbit-cli --test cli -- --list \
  | rg -Fqx "${rust_cli_permission_test}: test"; then
  printf 'error: release-soak CLI permission test is absent: %s\n' \
    "${rust_cli_permission_test}" >&2
  exit 1
fi
BSBIT_REQUIRE_PERMISSION_DENIAL=1 \
  /usr/bin/time -v -o "${stage_output}/rust-cli-process-time.txt" \
  cargo test --release --locked -p bsbit-cli --test cli \
  "${rust_cli_permission_test}" -- --exact \
  > "${stage_output}/rust-cli-process-result.txt" 2>&1

{
  printf 'schema\tbsbit-release-soak\n'
  printf 'captured_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'git_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'rustc\t%s\n' "$(rustc --version)"
  printf 'cargo\t%s\n' "$(cargo --version)"
  printf 'gcc\t%s\n' "$(gcc --version | head -n 1)"
  printf 'uname\t%s\n' "$(uname -a)"
  printf 'native_reader_cases\t%s\n' "${reader_cases}"
  printf 'native_record_cases\t%s\n' "${record_cases}"
  printf 'native_header_cases\t%s\n' "${header_cases}"
  printf 'rust_text_cases_per_test\t%s\n' "${text_cases}"
  printf 'cli_permission_denial_required\ttrue\n'
} > "${stage_output}/environment.txt"

sha256sum "${source_files[@]}" > "${stage_output}/source-inputs.post-run.sha256"
cmp -- "${stage_output}/source-inputs.pre-run.sha256" \
  "${stage_output}/source-inputs.post-run.sha256" \
  || die 'source inputs changed during release soak'
bash "${script_dir}/check-native-sources.sh" \
  > "${stage_output}/native-sources-after.txt"
cmp -- "${stage_output}/native-sources-before.txt" \
  "${stage_output}/native-sources-after.txt" \
  || die 'frozen sources changed during release soak'
(
  cd -- "${stage_output}"
  find . -type f ! -name SHA256SUMS -print0 | sort -z \
    | xargs -0 sha256sum > SHA256SUMS
)
mv -- "${stage_output}" "${final_output}"
completed=1
trap - EXIT
rm -rf -- "${native_root}"
printf 'Release soak captured at %s\n' "${final_output}"
