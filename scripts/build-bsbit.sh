#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

if (( $# > 1 )) || [[ ${1-} == --help || ${1-} == -h ]]; then
  printf 'usage: %s [OUTPUT_DIRECTORY]\n' "${0##*/}"
  exit $(( $# > 1 ? 2 : 0 ))
fi

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"
readonly output_root="$(realpath -m -- "${1:-${repository_root}/build/bsbit}")"
readonly target_dir="${output_root}/target"
readonly audit_dir="${output_root}/audit"
readonly binary="${output_root}/bsbit"
readonly pgo_input="${BSBIT_ALIGN_PGO_PROFILE-}"

for command in cargo cmp git install objdump python3 realpath rg rustc sha256sum uname; do
  command -v "${command}" >/dev/null 2>&1 || die "missing command: ${command}"
done
[[ -z ${RUSTFLAGS-} && -z ${CARGO_ENCODED_RUSTFLAGS-} ]] \
  || die 'unset RUSTFLAGS and CARGO_ENCODED_RUSTFLAGS'
[[ $(uname -m) == x86_64 ]] || die 'bsbit align requires x86-64'
rg -q -m 1 '(^|[[:space:]])sse4_2([[:space:]]|$)' /proc/cpuinfo \
  || die 'bsbit align requires SSE4.2'
rg -q -m 1 '(^|[[:space:]])avx2([[:space:]]|$)' /proc/cpuinfo \
  || die 'bsbit align requires AVX2'
rg -q -m 1 '(^|[[:space:]])popcnt([[:space:]]|$)' /proc/cpuinfo \
  || die 'bsbit align requires POPCNT'
case "${output_root}/" in
  "${repository_root}/build/"*|/tmp/*) ;;
  *) die 'output must be under repository build/ or /tmp' ;;
esac
[[ ! -e ${output_root} && ! -L ${output_root} ]] \
  || die "refusing to overwrite ${output_root}"
git -C "${repository_root}" diff --quiet --exit-code \
  || die 'tracked worktree changes exist; build audited artifacts only from committed source'
git -C "${repository_root}" diff --cached --quiet --exit-code \
  || die 'staged changes exist; build audited artifacts only from committed source'
[[ -z $(git -C "${repository_root}" ls-files --others --exclude-standard) ]] \
  || die 'untracked source files exist; build audited artifacts only from committed source'
readonly native_sources_before="$("${script_dir}/check-native-sources.sh")"

python3 "${repository_root}/scripts/check-release-notices.py" \
  --artifact binary --require-project-license

readonly base_rustflags='-C target-cpu=x86-64-v3 -C target-feature=+popcnt -A dead-code'
if [[ -n ${pgo_input} ]]; then
  readonly pgo_profile="$(realpath -e -- "${pgo_input}")"
  readonly rustflags="${base_rustflags} -C profile-use=${pgo_profile}"
else
  readonly pgo_profile=''
  readonly rustflags="${base_rustflags}"
fi

mkdir -p -- "$(dirname -- "${output_root}")"
mkdir -- "${output_root}" "${target_dir}" "${audit_dir}"
printf '%s\n' "${native_sources_before}" \
  > "${audit_dir}/native-sources-before.txt"
git -C "${repository_root}" status --short > "${audit_dir}/git-status-before-build.txt"
{
  printf 'captured_utc=%s\n' "$(date -u +%FT%TZ)"
  printf 'repository_commit=%s\n' "$(git -C "${repository_root}" rev-parse HEAD)"
  printf 'features=standard\n'
  printf 'rustflags=%s\n' "${rustflags}"
  printf 'pgo_profile=%s\n' "${pgo_profile:-none}"
  uname -a
  rustc --version --verbose
  cargo --version --verbose
} > "${audit_dir}/environment.txt"
if [[ -n ${pgo_profile} ]]; then
  sha256sum -- "${pgo_profile}" > "${audit_dir}/pgo-profile.sha256"
fi

(
  cd -- "${repository_root}"
  CARGO_INCREMENTAL=0 CARGO_TARGET_DIR="${target_dir}" RUSTFLAGS="${rustflags}" \
    CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
    /usr/bin/time -v cargo build --locked --release -p bsbit-cli --bin bsbit
) > "${audit_dir}/compile.stdout" 2> "${audit_dir}/compile.stderr"

install -m 0755 -- "${target_dir}/release/bsbit" "${binary}"
readonly popcnt_count="$(objdump -d -- "${binary}" | rg -c '[[:space:]]popcnt' || true)"
[[ ${popcnt_count} =~ ^[1-9][0-9]*$ ]] || die 'binary contains no POPCNT instruction'
readonly avx512_count="$(objdump -d -- "${binary}" | rg -c '(%zmm|%k[0-7])' || printf '0\n')"
readonly libsais_symbol_count="$(
  objdump -t -- "${binary}" \
    | rg -c '[[:space:]]libsais(64)?\.c$' \
    || true
)"
[[ ${libsais_symbol_count} =~ ^[1-9][0-9]*$ ]] \
  || die 'binary contains no linked libsais object symbol'
printf 'rust_isa_baseline=x86-64-v3\npopcnt_instruction_sites=%s\nwhole_elf_avx512_multiversion_sites=%s\nlibsais_object_symbols=%s\n' \
  "${popcnt_count}" "${avx512_count}" "${libsais_symbol_count}" \
  > "${audit_dir}/instruction-contract.txt"
sha256sum -- "${binary}" > "${audit_dir}/binary.sha256"
python3 "${repository_root}/scripts/check-release-notices.py" \
  --artifact binary --assemble "${output_root}/licenses" \
  > "${audit_dir}/license-assembly.txt"
[[ -f ${output_root}/licenses/third-party/libsais-2.10.4-LICENSE ]] \
  || die 'binary license assembly omitted libsais'
cmp -- \
  "${repository_root}/external/libsais/LICENSE" \
  "${output_root}/licenses/third-party/libsais-2.10.4-LICENSE" \
  || die 'assembled libsais license differs from the pinned source'
rg -q '"name"[[:space:]]*:[[:space:]]*"libsais"' \
  "${output_root}/licenses/license-manifest.json" \
  || die 'binary license manifest omitted libsais'
"${script_dir}/check-native-sources.sh" \
  > "${audit_dir}/native-sources-after.txt"
cmp -- \
  "${audit_dir}/native-sources-before.txt" \
  "${audit_dir}/native-sources-after.txt" \
  || die 'pinned native sources changed during build'
git -C "${repository_root}" status --short > "${audit_dir}/git-status-after-build.txt"
cmp -- "${audit_dir}/git-status-before-build.txt" "${audit_dir}/git-status-after-build.txt" \
  || die 'working tree changed during build'
printf 'audited bsbit binary built in %s\n' "${output_root}"
