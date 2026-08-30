#!/usr/bin/env bash
set -euo pipefail

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"

if [[ $# -ne 2 ]]; then
  printf 'usage: %s NORMAL-HTSLIB-PREFIX SANITIZED-HTSLIB-SOURCE\n' "$0" >&2
  exit 64
fi

readonly normal_prefix="$(realpath -- "$1")"
readonly sanitized_source="$(realpath -- "$2")"
readonly normal_pc="${normal_prefix}/lib/pkgconfig/htslib.pc"
readonly sanitized_archive="${sanitized_source}/libhts.a"

[[ -f ${normal_pc} ]] || {
  printf 'error: missing normal htslib.pc: %s\n' "${normal_pc}" >&2
  exit 1
}
[[ -f ${sanitized_source}/htslib/hts.h && -f ${sanitized_archive} ]] || {
  printf 'error: incomplete sanitized HTSlib source/build: %s\n' \
    "${sanitized_source}" >&2
  exit 1
}
[[ $(PKG_CONFIG_PATH="${normal_prefix}/lib/pkgconfig" pkg-config \
  --modversion htslib) == 1.24 ]] || {
  printf 'error: normal HTSlib pkg-config version is not exactly 1.24\n' >&2
  exit 1
}
rg -q '^#define HTS_VERSION 102400$' "${sanitized_source}/htslib/hts.h"
rg -q -- '-fsanitize=address,undefined' "${sanitized_source}/config.mk"

cd -- "${repository_root}"
readonly require_sites="$(rg -o 'REQUIRE\(' \
  crates/bsbit-hts/htslib-shim/tests/shim_smoke.c | wc -l)"
[[ ${require_sites} -eq 120 ]] || {
  printf 'error: expected 119 smoke assertions plus one macro definition, found %s sites\n' \
    "${require_sites}" >&2
  exit 1
}
readonly fault_require_sites="$(rg -o 'REQUIRE\(' \
  crates/bsbit-hts/htslib-shim/tests/fault_smoke.c | wc -l)"
[[ ${fault_require_sites} -eq 73 ]] || {
  printf 'error: expected 72 fault assertions plus one macro definition, found %s sites\n' \
    "${fault_require_sites}" >&2
  exit 1
}
readonly mutation_require_sites="$(rg -o 'REQUIRE\(' \
  crates/bsbit-hts/htslib-shim/tests/mutation_smoke.c | wc -l)"
[[ ${mutation_require_sites} -eq 42 ]] || {
  printf 'error: expected 41 mutation assertion sites plus one macro definition, found %s sites\n' \
    "${mutation_require_sites}" >&2
  exit 1
}
readonly scratch="$(mktemp -d --tmpdir bsbit-hts-shim.XXXXXX)"
readonly normal_build="${scratch}/normal"
readonly sanitizer_binary="${scratch}/sanitize/bsbit_htslib_shim_smoke"
readonly sanitizer_mutation_binary="${scratch}/sanitize/bsbit_htslib_shim_mutation_smoke"
readonly sanitizer_fault_binary="${scratch}/sanitize/bsbit_htslib_shim_fault_smoke"
cleanup() {
  rm -rf -- "${scratch}"
}
trap cleanup EXIT

PKG_CONFIG_PATH="${normal_prefix}/lib/pkgconfig" cmake \
  -S crates/bsbit-hts/htslib-shim -B "${normal_build}" -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "${normal_build}" --parallel 2
LD_LIBRARY_PATH="${normal_prefix}/lib" ctest \
  --test-dir "${normal_build}" --output-on-failure

mapfile -t actual_symbols < <(
  nm -g --defined-only "${normal_build}/libbsbit_htslib_shim.a" |
    awk 'NF == 3 { print $3 }' | sort -u
)
readonly expected_symbols=(
  bsbit_hts_bam_index_build
  bsbit_hts_bgzf_writer_destroy
  bsbit_hts_bgzf_writer_finish
  bsbit_hts_bgzf_writer_flush
  bsbit_hts_bgzf_writer_open
  bsbit_hts_bgzf_writer_write
  bsbit_hts_health_check
  bsbit_hts_indexed_reader_close
  bsbit_hts_indexed_reader_destroy
  bsbit_hts_indexed_reader_header_text
  bsbit_hts_indexed_reader_next
  bsbit_hts_indexed_reader_open
  bsbit_hts_indexed_reader_query
  bsbit_hts_indexed_reader_reference
  bsbit_hts_indexed_reader_reference_count
  bsbit_hts_indexed_fasta_reader_close
  bsbit_hts_indexed_fasta_reader_destroy
  bsbit_hts_indexed_fasta_reader_fetch
  bsbit_hts_indexed_fasta_reader_open
  bsbit_hts_indexed_fasta_reader_reference
  bsbit_hts_indexed_fasta_reader_reference_count
  bsbit_hts_pin_current_thread
  bsbit_hts_reader_close
  bsbit_hts_reader_compression
  bsbit_hts_reader_destroy
  bsbit_hts_reader_open
  bsbit_hts_reader_read
  bsbit_hts_runtime_version
  bsbit_hts_shim_abi_version
  bsbit_hts_tabix_index_build
  bsbit_hts_writer_destroy
  bsbit_hts_writer_finish
  bsbit_hts_writer_open_bam
  bsbit_hts_writer_open_bam_threads
  bsbit_hts_writer_open_bam_threads_level
  bsbit_hts_writer_write_bam_fields
  bsbit_hts_writer_write_record
)
mapfile -t sorted_expected < <(printf '%s\n' "${expected_symbols[@]}" | sort)
[[ ${actual_symbols[*]} == "${sorted_expected[*]}" ]] || {
  printf 'error: shim exported-symbol allow-list mismatch\n' >&2
  printf 'actual: %s\n' "${actual_symbols[*]}" >&2
  printf 'expected: %s\n' "${sorted_expected[*]}" >&2
  exit 1
}

mkdir -p -- "$(dirname -- "${sanitizer_binary}")"
gcc -std=c11 -D_XOPEN_SOURCE=700 \
  -Wall -Wextra -Wpedantic -Wconversion -Werror \
  -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -fno-sanitize-recover=undefined \
  -Icrates/bsbit-hts/htslib-shim -isystem "${sanitized_source}" \
  crates/bsbit-hts/htslib-shim/bsbit_hts.c \
  crates/bsbit-hts/htslib-shim/tests/shim_smoke.c \
  "${sanitized_archive}" -ldeflate -llzma -lbz2 -lz -lm -lpthread \
  -o "${sanitizer_binary}"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
  "${sanitizer_binary}" "${scratch}/sanitizer-artifacts"

gcc -std=c11 -D_XOPEN_SOURCE=700 \
  -Wall -Wextra -Wpedantic -Wconversion -Werror \
  -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -fno-sanitize-recover=undefined \
  -Icrates/bsbit-hts/htslib-shim -isystem "${sanitized_source}" \
  crates/bsbit-hts/htslib-shim/bsbit_hts.c \
  crates/bsbit-hts/htslib-shim/tests/mutation_smoke.c \
  "${sanitized_archive}" -ldeflate -llzma -lbz2 -lz -lm -lpthread \
  -o "${sanitizer_mutation_binary}"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
  "${sanitizer_mutation_binary}" "${scratch}/sanitizer-mutation-artifacts"

gcc -std=c11 -D_XOPEN_SOURCE=700 \
  -Wall -Wextra -Wpedantic -Wconversion -Werror \
  -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -fno-sanitize-recover=undefined \
  -Icrates/bsbit-hts/htslib-shim -isystem "${sanitized_source}" \
  crates/bsbit-hts/htslib-shim/bsbit_hts.c \
  crates/bsbit-hts/htslib-shim/tests/fault_smoke.c \
  "${sanitized_archive}" -ldeflate -llzma -lbz2 -lz -lm -lpthread \
  -Wl,--wrap=calloc \
  -Wl,--wrap=malloc \
  -Wl,--wrap=bgzf_open \
  -Wl,--wrap=bgzf_compression \
  -Wl,--wrap=bgzf_read \
  -Wl,--wrap=bgzf_write \
  -Wl,--wrap=bgzf_flush \
  -Wl,--wrap=bgzf_close \
  -Wl,--wrap=sam_hdr_parse \
  -Wl,--wrap=bam_init1 \
  -Wl,--wrap=hts_open \
  -Wl,--wrap=sam_hdr_write \
  -Wl,--wrap=sam_parse1 \
  -Wl,--wrap=sam_write1 \
  -Wl,--wrap=hts_close \
  -o "${sanitizer_fault_binary}"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
  "${sanitizer_fault_binary}" "${scratch}/sanitizer-fault-artifacts"

"${script_dir}/check-native-sources.sh"
printf 'HTSlib shim checks passed\n'
