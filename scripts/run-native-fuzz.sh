#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

readonly expected_clang='Ubuntu clang version 18.1.3 (1ubuntu1)'
readonly expected_llvm='18.1.3'
readonly fuzz_seed=42574257
readonly max_len=4096
readonly max_decoded_bytes=1048576
readonly targets=(reader record header)
readonly requested_output=${1:-}
readonly seconds_per_target=${BSBIT_NATIVE_FUZZ_SECONDS_PER_TARGET:-30}
readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(git -C "${script_dir}/.." rev-parse --show-toplevel)"

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

if (( $# > 1 )); then
  printf 'usage: %s [ABSENT_OUTPUT_DIRECTORY]\n' "$0" >&2
  exit 2
fi
if [[ ! ${seconds_per_target} =~ ^[0-9]+$ ]] \
  || (( seconds_per_target < 1 || seconds_per_target > 3600 )); then
  die 'BSBIT_NATIVE_FUZZ_SECONDS_PER_TARGET must be an integer from 1 through 3600'
fi

cd -- "${repository_root}"
for command in ar autoreconf awk clang cmp cp find git gzip llvm-config make \
  mkdir mktemp mv nm python3 rg rm sed sha256sum sort tr uname wc xargs; do
  command -v "${command}" >/dev/null 2>&1 \
    || die "required native fuzz tool is unavailable: ${command}"
done
[[ $(clang --version | sed -n '1p') == "${expected_clang}" ]] \
  || die "clang must be exactly ${expected_clang}"
[[ $(llvm-config --version) == "${expected_llvm}" ]] \
  || die "LLVM tools must be exactly ${expected_llvm}"

formal_capture=0
ephemeral_output=0
if [[ -n ${requested_output} ]]; then
  formal_capture=1
  [[ ! -e ${requested_output} && ! -L ${requested_output} ]] \
    || die "output already exists: ${requested_output}"
  output_parent="$(dirname -- "${requested_output}")"
  [[ -d ${output_parent} ]] || die "output parent does not exist: ${output_parent}"
  output_parent="$(cd -- "${output_parent}" && pwd -P)"
  final_output="${output_parent}/$(basename -- "${requested_output}")"
  git diff --quiet --exit-code \
    || die 'tracked worktree changes exist; capture only from committed source'
  git diff --cached --quiet --exit-code \
    || die 'staged changes exist; capture only from committed source'
  [[ -z $(git ls-files --others --exclude-standard) ]] \
    || die 'untracked source files exist; capture only from committed source'
  for tracked_input in \
    crates/bsbit-hts/htslib-shim/bsbit_hts.h \
    crates/bsbit-hts/htslib-shim/bsbit_hts.c \
    crates/bsbit-hts/htslib-shim/tests/fuzz/README.md \
    crates/bsbit-hts/htslib-shim/tests/fuzz/native_fuzz.c \
    scripts/run-native-fuzz.sh \
    scripts/summarize-native-fuzz.py
  do
    git ls-files --error-unmatch -- "${tracked_input}" >/dev/null
  done
  stage_output="$(mktemp -d -- "${output_parent}/.$(basename -- "${requested_output}").partial.XXXXXX")"
else
  ephemeral_output=1
  final_output=''
  stage_output="$(mktemp -d --tmpdir bsbit-native-fuzz-output.XXXXXX)"
fi
native_root="$(mktemp -d --tmpdir bsbit-native-fuzz-build.XXXXXX)"
native_tmp="${native_root}/tmp"
mkdir -- "${native_tmp}"
readonly formal_capture ephemeral_output final_output stage_output native_root native_tmp

remove_owned_tree() {
  local path=$1
  local prefix=$2
  [[ -d ${path} ]] || return 0
  if [[ ${path##*/} != "${prefix}."* ]]; then
    printf 'refusing to remove unexpected temporary tree: %s\n' "${path}" >&2
    return 0
  fi
  rm -rf -- "${path}"
}

cleanup() {
  local status=$?
  if (( status != 0 )); then
    printf '%s\n' "${status}" > "${stage_output}/exit-status.txt"
    if [[ -f ${native_root}/htslib/config.log ]]; then
      cp -- "${native_root}/htslib/config.log" \
        "${stage_output}/htslib-config.log" || true
    fi
    scripts/check-native-sources.sh \
      > "${stage_output}/external-after-failure.txt" 2>&1 || true
    printf 'native fuzz failure evidence retained at %s\n' "${stage_output}" >&2
  elif (( ephemeral_output == 1 )); then
    remove_owned_tree "${stage_output}" bsbit-native-fuzz-output
  fi
  remove_owned_tree "${native_root}" bsbit-native-fuzz-build
  trap - EXIT
  exit "${status}"
}
trap cleanup EXIT

readonly -a source_files=(
  .gitmodules
  crates/bsbit-hts/htslib-shim/bsbit_hts.h
  crates/bsbit-hts/htslib-shim/bsbit_hts.c
  crates/bsbit-hts/htslib-shim/tests/fuzz/README.md
  crates/bsbit-hts/htslib-shim/tests/fuzz/native_fuzz.c
  scripts/run-native-fuzz.sh
  scripts/check-native-sources.sh
  scripts/summarize-native-fuzz.py
)

mkdir -p \
  "${stage_output}/artifacts/reader" \
  "${stage_output}/artifacts/record" \
  "${stage_output}/artifacts/header" \
  "${stage_output}/corpus/reader" \
  "${stage_output}/corpus/record" \
  "${stage_output}/corpus/header" \
  "${stage_output}/logs" \
  "${native_root}/objects" \
  "${native_root}/scratch/reader" \
  "${native_root}/scratch/record" \
  "${native_root}/scratch/header" \
  "${native_root}/htslib"

sha256sum "${source_files[@]}" > "${stage_output}/source-before.sha256"
scripts/check-native-sources.sh > "${stage_output}/native-sources-before.txt"
git rev-parse HEAD > "${stage_output}/commit.txt"
{
  printf 'schema=bsbit-native-coverage-fuzz\n'
  printf 'clang=%s\n' "$(clang --version | tr '\n' ';')"
  printf 'llvm=%s\n' "$(llvm-config --version)"
  printf 'autoreconf=%s\n' "$(autoreconf --version | sed -n '1p')"
  printf 'make=%s\n' "$(make --version | sed -n '1p')"
  printf 'gzip=%s\n' "$(gzip --version | sed -n '1p')"
  printf 'python=%s\n' "$(python3 --version)"
  printf 'kernel=%s\n' "$(uname -a)"
  printf 'seconds_per_target=%s\n' "${seconds_per_target}"
  printf 'seed=%s\n' "${fuzz_seed}"
  printf 'max_len=%s\n' "${max_len}"
  printf 'max_decoded_bytes=%s\n' "${max_decoded_bytes}"
  printf 'asan_options=detect_leaks=0:halt_on_error=1:abort_on_error=1:symbolize=1\n'
  printf 'ubsan_options=halt_on_error=1:print_stacktrace=1\n'
  printf 'htslib_commit=%s\n' "$(git -C external/htslib rev-parse HEAD)"
  printf 'htscodecs_commit=%s\n' \
    "$(git -C external/htslib/htscodecs rev-parse HEAD)"
  printf 'libsais_commit=%s\n' "$(git -C external/libsais rev-parse HEAD)"
} > "${stage_output}/environment.txt"

cp -a -- external/htslib/. "${native_root}/htslib/"
rm -rf -- "${native_root}/htslib/.git" \
  "${native_root}/htslib/htscodecs/.git"
printf '#define HTSCODECS_VERSION_TEXT "1.6.7"\n' \
  > "${native_root}/htslib/htscodecs/htscodecs/version.h"

export TMPDIR="${native_tmp}"
export TMP="${native_tmp}"
export TEMP="${native_tmp}"
export ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:abort_on_error=1:symbolize=1'
export UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'
(
  cd -- "${native_root}/htslib"
  autoreconf -i
  CC=clang \
  CFLAGS='-O1 -g -fno-omit-frame-pointer -fno-sanitize-recover=undefined -fsanitize=fuzzer-no-link,address,undefined' \
  LDFLAGS='-fsanitize=address,undefined' \
    ./configure --disable-plugins --disable-libcurl --disable-gcs --disable-s3
  make -j2 lib-static
) > "${stage_output}/htslib-build.log" 2>&1

readonly -a compile_flags=(
  -std=c11
  -D_XOPEN_SOURCE=700
  -Wall
  -Wextra
  -Wpedantic
  -Wconversion
  -Werror
  -O1
  -g
  -fno-omit-frame-pointer
  -fno-sanitize-recover=undefined
  -fsanitize=fuzzer-no-link,address,undefined
  -Icrates/bsbit-hts/htslib-shim
  -isystem
  "${native_root}/htslib"
)
clang "${compile_flags[@]}" -c crates/bsbit-hts/htslib-shim/bsbit_hts.c \
  -o "${native_root}/objects/bsbit_hts.o"
for target_spec in reader:1 record:2 header:3; do
  target_name=${target_spec%%:*}
  target_mode=${target_spec##*:}
  clang "${compile_flags[@]}" -DBSBIT_NATIVE_FUZZ_MODE="${target_mode}" \
    -c crates/bsbit-hts/htslib-shim/tests/fuzz/native_fuzz.c \
    -o "${native_root}/objects/${target_name}.o"
  clang -O1 -g -fno-omit-frame-pointer -fno-sanitize-recover=undefined \
    -fsanitize=fuzzer,address,undefined \
    "${native_root}/objects/bsbit_hts.o" \
    "${native_root}/objects/${target_name}.o" \
    "${native_root}/htslib/libhts.a" \
    -ldeflate -llzma -lbz2 -lz -lm -lpthread \
    -o "${native_root}/${target_name}_fuzz"
done
sha256sum "${native_root}/reader_fuzz" "${native_root}/record_fuzz" \
  "${native_root}/header_fuzz" > "${stage_output}/binaries.sha256"

ar t "${native_root}/htslib/libhts.a" > "${stage_output}/archive-members.txt"
nm -A --undefined-only "${native_root}/htslib/libhts.a" \
  > "${stage_output}/instrumentation-symbols.txt"
archive_members=$(wc -l < "${stage_output}/archive-members.txt")
coverage_members=$(awk '/__sanitizer_cov_8bit_counters_init/ {print $1}' \
  "${stage_output}/instrumentation-symbols.txt" | sort -u | wc -l)
asan_members=$(awk '/__asan_init/ {print $1}' \
  "${stage_output}/instrumentation-symbols.txt" | sort -u | wc -l)
ubsan_members=$(awk '/__ubsan_handle_/ {print $1}' \
  "${stage_output}/instrumentation-symbols.txt" | sort -u | wc -l)
[[ ${archive_members} -eq 52 ]]
[[ ${coverage_members} -eq 51 ]]
[[ ${asan_members} -eq 52 ]]
[[ ${ubsan_members} -eq 50 ]]
rg -q 'libhts\.a:bgzf\.o:.*__sanitizer_cov_8bit_counters_init' \
  "${stage_output}/instrumentation-symbols.txt"
rg -q 'libhts\.a:htscodecs\.o:.*__sanitizer_cov_8bit_counters_init' \
  "${stage_output}/instrumentation-symbols.txt"
for object in bsbit_hts reader record header; do
  nm --undefined-only "${native_root}/objects/${object}.o" \
    > "${stage_output}/${object}-symbols.txt"
  rg -q '__sanitizer_cov_8bit_counters_init' \
    "${stage_output}/${object}-symbols.txt"
  rg -q '__asan_init' "${stage_output}/${object}-symbols.txt"
  rg -q '__ubsan_handle_' "${stage_output}/${object}-symbols.txt"
done
{
  printf 'component\ttotal_objects\tcoverage_objects\tasan_objects\tubsan_objects\n'
  printf 'htslib-static\t%s\t%s\t%s\t%s\n' \
    "${archive_members}" "${coverage_members}" "${asan_members}" "${ubsan_members}"
  printf 'project-shim\t1\t1\t1\t1\n'
  printf 'reader-target\t1\t1\t1\t1\n'
  printf 'record-target\t1\t1\t1\t1\n'
  printf 'header-target\t1\t1\t1\t1\n'
} > "${stage_output}/instrumentation.tsv"

printf '@r\nACGT\n+\nABCD\n' > "${stage_output}/corpus/reader/plain"
gzip -n -c "${stage_output}/corpus/reader/plain" \
  > "${stage_output}/corpus/reader/gzip"
printf '\037\213\010\000truncated' > "${stage_output}/corpus/reader/truncated-gzip"
printf 'r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tABCD\n' \
  > "${stage_output}/corpus/record/valid"
printf 'not-a-sam-record\n' > "${stage_output}/corpus/record/invalid"
printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n' \
  > "${stage_output}/corpus/header/valid"
printf '@SQ\tSN:missing-length\n' > "${stage_output}/corpus/header/invalid"

python3 scripts/summarize-native-fuzz.py self-test
for target_name in "${targets[@]}"; do
  BSBIT_NATIVE_FUZZ_SCRATCH="${native_root}/scratch/${target_name}" \
    "${native_root}/${target_name}_fuzz" \
    "${stage_output}/corpus/${target_name}" \
    "-max_total_time=${seconds_per_target}" \
    "-max_len=${max_len}" \
    -timeout=5 \
    -rss_limit_mb=2048 \
    "-seed=${fuzz_seed}" \
    -print_final_stats=1 \
    "-artifact_prefix=${stage_output}/artifacts/${target_name}/" \
    > "${stage_output}/logs/${target_name}.log" 2>&1
  if [[ -n $(find "${stage_output}/artifacts/${target_name}" -type f -print -quit) ]]; then
    die "native fuzz target produced a failure artifact: ${target_name}"
  fi
done

python3 scripts/summarize-native-fuzz.py summary \
  --requested-seconds "${seconds_per_target}" \
  --seed "${fuzz_seed}" \
  --max-len "${max_len}" \
  --output "${stage_output}/summary.tsv" \
  "${stage_output}/logs/reader.log" \
  "${stage_output}/logs/record.log" \
  "${stage_output}/logs/header.log"

sha256sum "${source_files[@]}" > "${stage_output}/source-after.sha256"
cmp -- "${stage_output}/source-before.sha256" "${stage_output}/source-after.sha256"
scripts/check-native-sources.sh > "${stage_output}/native-sources-after.txt"
cmp -- "${stage_output}/native-sources-before.txt" \
  "${stage_output}/native-sources-after.txt"
printf 'PASS\n' > "${stage_output}/result.txt"
(
  cd -- "${stage_output}"
  find . -type f ! -name SHA256SUMS -printf '%P\0' \
    | sort -z \
    | xargs -0 sha256sum \
    > SHA256SUMS
)

if (( formal_capture == 1 )); then
  mv -- "${stage_output}" "${final_output}"
  printf 'Native coverage fuzz artifact published at %s\n' "${final_output}"
  sed -n '1,4p' "${final_output}/summary.tsv"
else
  sed -n '1,4p' "${stage_output}/summary.tsv"
fi
