#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage: tests/qualification/check-simd-backends.sh [BSBIT_BINARY]

Force index construction and all four alignment modes (single/paired times
default/sensitive) through every backend supported by the current CPU. Require
byte-identical BAM output, and prove that unavailable backends fail before
input processing and leave empty direct outputs. Temporary evidence is removed
only on success.
EOF
}

if [[ ${1-} == --help || ${1-} == -h ]]; then
  usage
  exit 0
fi
[[ $# -le 1 ]] || die 'expected at most one binary path'

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/../.." && pwd)"
readonly binary="$(realpath -e -- "${1:-${repository_root}/target/release/bsbit}")"
[[ -x ${binary} ]] || die "binary is not executable: ${binary}"
for command in cmp grep mktemp realpath samtools sed; do
  command -v "${command}" >/dev/null 2>&1 || die "missing command: ${command}"
done

scratch="$(mktemp -d --tmpdir bsbit-simd-backends.XXXXXX)"
readonly scratch
completed=0
cleanup() {
  if [[ ${completed} == 1 ]]; then
    rm -rf -- "${scratch}"
  else
    printf 'SIMD differential evidence retained at %s\n' "${scratch}" >&2
  fi
}
trap cleanup EXIT

cp -- "${repository_root}/docs/examples/test.fa" "${scratch}/reference.fa"
samtools faidx "${scratch}/reference.fa"

"${binary}" cpu > "${scratch}/cpu.txt"
selected="$(sed -n 's/^backend=//p' "${scratch}/cpu.txt")"
[[ -n ${selected} ]] || die 'CPU diagnostic did not report a backend'
architecture="$(sed -n 's/^architecture=//p' "${scratch}/cpu.txt")"
[[ -n ${architecture} ]] || die 'CPU diagnostic did not report an architecture'

case "${architecture}" in
  x86_64)
    available_backends=(scalar)
    unavailable_backends=(neon)
    if grep -qx 'sse2=1' "${scratch}/cpu.txt"; then
      available_backends+=(sse2)
    else
      unavailable_backends+=(sse2)
    fi
    if grep -qx 'sse4.1=1' "${scratch}/cpu.txt" \
      && grep -qx 'sse4.2=1' "${scratch}/cpu.txt" \
      && grep -qx 'popcnt=1' "${scratch}/cpu.txt"; then
      available_backends+=(sse4.2)
    else
      unavailable_backends+=(sse4.2)
    fi
    if grep -qx 'avx2=1' "${scratch}/cpu.txt" \
      && grep -qx 'popcnt=1' "${scratch}/cpu.txt"; then
      available_backends+=(avx2)
    else
      unavailable_backends+=(avx2)
    fi
    if grep -qx 'avx512f=1' "${scratch}/cpu.txt" \
      && grep -qx 'avx512bw=1' "${scratch}/cpu.txt" \
      && grep -qx 'avx2=1' "${scratch}/cpu.txt" \
      && grep -qx 'popcnt=1' "${scratch}/cpu.txt"; then
      available_backends+=(avx512)
    else
      unavailable_backends+=(avx512)
    fi
    ;;
  aarch64)
    available_backends=(scalar)
    unavailable_backends=(sse2 sse4.2 avx2 avx512)
    if grep -qx 'neon=1' "${scratch}/cpu.txt"; then
      available_backends+=(neon)
    else
      unavailable_backends+=(neon)
    fi
    ;;
  *)
    available_backends=(scalar)
    unavailable_backends=(sse2 sse4.2 avx2 avx512 neon)
    ;;
esac

build_index() {
  local backend=$1
  local suffix=$2
  "${binary}" index \
    --reference "${scratch}/reference.fa" \
    --output "${scratch}/reference-${suffix}.bsbit" \
    --threads 1 \
    --memory-mib 512 \
    --simd-backend "${backend}"
}

align_mode() {
  local backend=$1
  local suffix=$2
  local mode=$3
  local index_path="${scratch}/reference-${suffix}.bsbit"
  local output_path="${scratch}/${mode}-${suffix}.bam"
  local arguments=(
    align
    --index "${index_path}"
    --read1 "${repository_root}/docs/examples/test_R1.fastq"
    --output "${output_path}"
    --threads 1
    --compression-threads 0
    --compression-level 0
    --simd-backend "${backend}"
  )
  case "${mode}" in
    single-default) ;;
    single-sensitive) arguments+=(--sensitive) ;;
    paired-default)
      arguments+=(
        --read2 "${repository_root}/docs/examples/test_R2.fastq"
        --min-template-span 100
        --max-template-span 250
      )
      ;;
    paired-sensitive)
      arguments+=(
        --read2 "${repository_root}/docs/examples/test_R2.fastq"
        --min-template-span 100
        --max-template-span 250
        --sensitive
      )
      ;;
    *) die "internal unsupported alignment mode: ${mode}" ;;
  esac
  "${binary}" "${arguments[@]}"
  samtools quickcheck -v "${output_path}"
}

build_index auto auto

for backend in "${unavailable_backends[@]}"; do
  if "${binary}" cpu --simd-backend "${backend}" \
    > "${scratch}/cpu-${backend}.stdout" \
    2> "${scratch}/cpu-${backend}.stderr"; then
    die "unsupported ${backend} backend was accepted"
  fi
  grep -q 'requires ' "${scratch}/cpu-${backend}.stderr" \
    || die "unsupported ${backend} CPU error did not explain its requirement"

  if "${binary}" index \
    --reference "${scratch}/reference.fa" \
    --output "${scratch}/unsupported-index-${backend}.bsbit" \
    --threads 1 \
    --memory-mib 512 \
    --simd-backend "${backend}" \
    > "${scratch}/index-${backend}.stdout" \
    2> "${scratch}/index-${backend}.stderr"; then
    die "index construction accepted unsupported ${backend} backend"
  fi
  grep -q 'requires ' "${scratch}/index-${backend}.stderr" \
    || die "index's unsupported ${backend} error did not explain its requirement"
  [[ -f "${scratch}/unsupported-index-${backend}.bsbit" \
    && ! -s "${scratch}/unsupported-index-${backend}.bsbit" ]] \
    || die "unsupported ${backend} index construction did not leave an empty output"

  if "${binary}" align \
    --index "${scratch}/reference-auto.bsbit" \
    --read1 "${repository_root}/docs/examples/test_R1.fastq" \
    --output "${scratch}/unsupported-align-${backend}.bam" \
    --simd-backend "${backend}" \
    > "${scratch}/align-${backend}.stdout" \
    2> "${scratch}/align-${backend}.stderr"; then
    die "alignment accepted unsupported ${backend} backend"
  fi
  grep -q 'requires ' "${scratch}/align-${backend}.stderr" \
    || die "alignment's unsupported ${backend} error did not explain its requirement"
  [[ -f "${scratch}/unsupported-align-${backend}.bam" \
    && ! -s "${scratch}/unsupported-align-${backend}.bam" ]] \
    || die "unsupported ${backend} alignment did not leave an empty output"
done

readonly modes=(single-default single-sensitive paired-default paired-sensitive)
for mode in "${modes[@]}"; do
  align_mode auto auto "${mode}"
done

for backend in "${available_backends[@]}"; do
  "${binary}" cpu --simd-backend "${backend}" > "${scratch}/cpu-${backend}.txt"
  build_index "${backend}" "${backend}"
  for mode in "${modes[@]}"; do
    align_mode "${backend}" "${backend}" "${mode}"
  done
done

readonly baseline_backend="${available_backends[0]}"
for mode in "${modes[@]}"; do
  expected="${scratch}/${mode}-${baseline_backend}.bam"
  for backend in "${available_backends[@]}"; do
    cmp -- "${expected}" "${scratch}/${mode}-${backend}.bam" \
      || die "${backend} ${mode} BAM differs from ${baseline_backend}"
  done
  cmp -- "${expected}" "${scratch}/${mode}-auto.bam" \
    || die "automatic ${selected} ${mode} BAM differs from ${baseline_backend}"
done

printf 'SIMD backend differential passed (single/paired x default/sensitive): auto=%s; forced=%s\n' \
  "${selected}" "${available_backends[*]}"
completed=1
