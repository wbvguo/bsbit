#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

usage() {
  printf 'usage: %s [--verify-current] [OUTPUT_DIRECTORY]\n' "${0##*/}"
}

verify_current=0
if [[ ${1-} == --verify-current ]]; then
  verify_current=1
  shift
fi
if [[ ${1-} == --help || ${1-} == -h ]]; then
  usage
  exit 0
fi
if (( $# > 1 )); then
  usage >&2
  exit 2
fi
[[ ${1-} != -* ]] || die "unknown option: ${1}"
readonly verify_current

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"
if (( verify_current )); then
  readonly build_mode=current-working-tree-verification
  readonly release_eligible=no
  readonly default_output="${repository_root}/build/bsbit-verification"
  readonly binary_name=bsbit-NOT-FOR-RELEASE
else
  readonly build_mode=release
  readonly release_eligible=yes
  readonly default_output="${repository_root}/build/bsbit"
  readonly binary_name=bsbit
fi
readonly output_root="$(realpath -m -- "${1:-${default_output}}")"
readonly target_dir="${output_root}/target"
readonly audit_dir="${output_root}/audit"
readonly binary="${output_root}/${binary_name}"
readonly pgo_input="${BSBIT_ALIGN_PGO_PROFILE-}"

for command in awk cargo cmp git install objdump python3 realpath rg rustc sed sha256sum sort uname; do
  command -v "${command}" >/dev/null 2>&1 || die "missing command: ${command}"
done
for variable in \
  RUSTFLAGS CARGO_ENCODED_RUSTFLAGS \
  CARGO_BUILD_TARGET CARGO_BUILD_RUSTFLAGS \
  CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS \
  CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_RUSTFLAGS \
  CFLAGS CXXFLAGS CPPFLAGS \
  TARGET_CFLAGS TARGET_CXXFLAGS TARGET_CPPFLAGS \
  HOST_CFLAGS HOST_CXXFLAGS HOST_CPPFLAGS \
  CFLAGS_x86_64_unknown_linux_gnu \
  CXXFLAGS_x86_64_unknown_linux_gnu \
  CPPFLAGS_x86_64_unknown_linux_gnu \
  CFLAGS_aarch64_unknown_linux_gnu \
  CXXFLAGS_aarch64_unknown_linux_gnu \
  CPPFLAGS_aarch64_unknown_linux_gnu; do
  [[ -z ${!variable-} ]] || die "unset ${variable}"
done
readonly machine="$(uname -m)"
case "${machine}" in
  x86_64) readonly architecture=x86_64 ;;
  aarch64|arm64) readonly architecture=aarch64 ;;
  *) die "unsupported build architecture: ${machine}" ;;
esac
case "${output_root}/" in
  "${repository_root}/build/"*|/tmp/*) ;;
  *) die 'output must be under repository build/ or /tmp' ;;
esac
[[ ! -e ${output_root} && ! -L ${output_root} ]] \
  || die "refusing to overwrite ${output_root}"
if (( verify_current )); then
  printf '%s\n' \
    'warning: --verify-current builds a non-release artifact from the current working tree' \
    >&2
  if repository_commit="$(git -C "${repository_root}" rev-parse --verify HEAD 2>/dev/null)"; then
    git_metadata=available
  else
    repository_commit=unavailable
    git_metadata=unavailable
    printf '%s\n' \
      'warning: Git worktree metadata is unavailable; source identity cannot be audited' \
      >&2
  fi
  if native_sources_before="$("${script_dir}/check-native-sources.sh" 2>&1)"; then
    native_source_identity=verified
  else
    native_source_identity=unverified
    native_sources_before="unavailable: check-native-sources.sh failed
${native_sources_before}"
    printf '%s\n' \
      'warning: pinned native-source identity is unavailable in this working tree' \
      >&2
  fi
else
  git -C "${repository_root}" rev-parse --verify HEAD >/dev/null 2>&1 \
    || die 'Git worktree metadata is unavailable; release source cannot be audited'
  git -C "${repository_root}" diff --quiet --exit-code \
    || die 'tracked worktree changes exist; build audited artifacts only from committed source'
  git -C "${repository_root}" diff --cached --quiet --exit-code \
    || die 'staged changes exist; build audited artifacts only from committed source'
  [[ -z $(git -C "${repository_root}" ls-files --others --exclude-standard) ]] \
    || die 'untracked source files exist; build audited artifacts only from committed source'
  repository_commit="$(git -C "${repository_root}" rev-parse HEAD)"
  git_metadata=available
  native_sources_before="$("${script_dir}/check-native-sources.sh")"
  native_source_identity=verified
fi
readonly repository_commit git_metadata native_sources_before native_source_identity

bash "${script_dir}/check-release-source.sh"

python3 "${repository_root}/scripts/check-release-notices.py" \
  --artifact binary --require-project-license

# Keep ordinary Rust functions at the supported architecture's generic
# baseline. x86 extensions remain in guarded leaves; standard AArch64 Linux
# includes NEON in its architectural Rust target baseline.
case "${architecture}" in
  x86_64)
    readonly target_cpu=x86-64
    readonly expected_baseline_target_features=$'target_feature="fxsr"\ntarget_feature="sse"\ntarget_feature="sse2"'
    ;;
  aarch64)
    readonly target_cpu=generic
    readonly expected_baseline_target_features='target_feature="neon"'
    ;;
esac
readonly base_rustflags="-C target-cpu=${target_cpu}"
readonly baseline_target_features="$(rustc --print cfg -C target-cpu="${target_cpu}" \
  | rg '^target_feature=' | sort)"
[[ ${baseline_target_features} == "${expected_baseline_target_features}" ]] \
  || die "build baseline exposes unexpected ${architecture} target features"
if [[ -n ${pgo_input} ]]; then
  readonly pgo_profile="$(realpath -e -- "${pgo_input}")"
  readonly rustflags="${base_rustflags} -C profile-use=${pgo_profile}"
else
  readonly pgo_profile=''
  readonly rustflags="${base_rustflags}"
fi

mkdir -p -- "$(dirname -- "${output_root}")"
mkdir -- "${output_root}" "${target_dir}" "${audit_dir}"
if (( verify_current )); then
  printf '%s\n' \
    'This artifact was built with --verify-current and is NOT FOR RELEASE.' \
    > "${output_root}/NOT-FOR-RELEASE"
fi
printf 'build_mode=%s\nrelease_eligible=%s\ngit_metadata=%s\nnative_source_identity=%s\n' \
  "${build_mode}" "${release_eligible}" "${git_metadata}" \
  "${native_source_identity}" > "${audit_dir}/artifact-status.txt"
printf '%s\n' "${baseline_target_features}" \
  > "${audit_dir}/baseline-target-features.txt"
printf '%s\n' "${native_sources_before}" \
  > "${audit_dir}/native-sources-before.txt"
if [[ ${git_metadata} == available ]]; then
  git -C "${repository_root}" status --short \
    > "${audit_dir}/git-status-before-build.txt"
else
  printf '%s\n' 'unavailable: Git worktree metadata could not be resolved' \
    > "${audit_dir}/git-status-before-build.txt"
fi
{
  printf 'captured_utc=%s\n' "$(date -u +%FT%TZ)"
  printf 'build_mode=%s\n' "${build_mode}"
  printf 'release_eligible=%s\n' "${release_eligible}"
  printf 'repository_commit=%s\n' "${repository_commit}"
  printf 'git_metadata=%s\n' "${git_metadata}"
  printf 'native_source_identity=%s\n' "${native_source_identity}"
  printf 'architecture=%s\n' "${architecture}"
  printf 'target_cpu=%s\n' "${target_cpu}"
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
if [[ ${architecture} == x86_64 ]]; then
readonly popcnt_count="$(objdump -d -- "${binary}" | rg -c '[[:space:]]popcnt' || true)"
[[ ${popcnt_count} =~ ^[1-9][0-9]*$ ]] || die 'binary contains no POPCNT instruction'
readonly avx_count="$(objdump -d -- "${binary}" | rg -c '(%ymm|[[:space:]]v[a-z0-9]+)' || true)"
[[ ${avx_count} =~ ^[1-9][0-9]*$ ]] || die 'binary contains no AVX2 backend instruction'
readonly avx512_count="$(objdump -d -- "${binary}" | rg -c '(%zmm|%k[0-7])' || printf '0\n')"
[[ ${avx512_count} =~ ^[1-9][0-9]*$ ]] || die 'binary contains no AVX-512 backend instruction'
readonly sse_above_baseline_pattern='[[:space:]](addsubp[ds]|haddp[ds]|hsubp[ds]|lddqu|movddup|movshdup|movsldup|fisttp[lqs]?|monitor|mwait|pabs[bwd]|palignr|phadd[wd]|phaddsw|phsub[wd]|phsubsw|pmaddubsw|pmulhrsw|pshufb|psign[bwd]|blendp[ds]|blendvp[ds]|dpp[ds]|extractps|insertps|movntdqa|mpsadbw|packusdw|pblendvb|pblendw|pcmpeqq|pextr[bdq]|phminposuw|pinsr[bdq]|pmaxs[bd]|pmaxu[dw]|pmins[bd]|pminu[dw]|pmovsx(bd|bq|bw|dq|dw|wq)|pmovzx(bd|bq|bw|dq|dw|wq)|pmuldq|pmulld|ptest|roundp[ds]|rounds[ds]|crc32|pcmpestr[im]|pcmpistr[im]|pcmpgtq)([[:space:]]|$)'
readonly sse_above_baseline_count="$(
  objdump -d -- "${binary}" | rg -c "${sse_above_baseline_pattern}" || printf '0\n'
)"
# AVX2 does not imply ABM/BMI, FMA, AES/PCLMUL, F16C, ADX, MOVBE, or
# random-number instructions. None belongs to bsbit's guarded backend contract.
readonly uncontracted_instruction_pattern='[[:space:]](lzcnt|andn|bextr|blsi|blsmsk|blsr|bzhi|pdep|pext|mulx|rorx|sarx|shlx|shrx|vfmadd[^[:space:]]*|vfmsub[^[:space:]]*|vfnmadd[^[:space:]]*|vfnmsub[^[:space:]]*|v?aesenc|v?aesenclast|v?aesdec|v?aesdeclast|v?aesimc|v?aeskeygenassist|v?pclmulqdq|v?gf2p8[^[:space:]]*|vcvtph2ps|vcvtps2ph|adcx|adox|movbe|rdrand|rdseed)([[:space:]]|$)'
readonly uncontracted_instruction_count="$(
  objdump -d -- "${binary}" | rg -c "${uncontracted_instruction_pattern}" || printf '0\n'
)"
[[ ${uncontracted_instruction_count} == 0 ]] \
  || die 'binary contains an instruction outside the x86-64 SIMD contract'

readonly popcnt_owners="${audit_dir}/popcnt-instruction-owners.txt"
readonly avx_owners="${audit_dir}/avx-instruction-owners.txt"
readonly avx512_owners="${audit_dir}/avx512-instruction-owners.txt"
readonly sse_above_baseline_owners="${audit_dir}/sse-above-baseline-instruction-owners.txt"
readonly sha_owners="${audit_dir}/sha-instruction-owners.txt"
readonly xgetbv_owners="${audit_dir}/xgetbv-instruction-owners.txt"
objdump -dC -- "${binary}" \
  | sed -n '/^[[:xdigit:]][[:xdigit:]]* <.*>:/h; /[[:space:]]popcnt/{x;p;x;}' \
  | sort -u > "${popcnt_owners}"
objdump -dC -- "${binary}" \
  | sed -n '/^[[:xdigit:]][[:xdigit:]]* <.*>:/h; /%ymm\|[[:space:]]v[a-z0-9][a-z0-9]*/{x;p;x;}' \
  | sort -u > "${avx_owners}"
objdump -dC -- "${binary}" \
  | sed -n '/^[[:xdigit:]][[:xdigit:]]* <.*>:/h; /%zmm\|%k[0-7]/{x;p;x;}' \
  | sort -u > "${avx512_owners}"
objdump -dC -- "${binary}" \
  | awk -v pattern="${sse_above_baseline_pattern}" \
      '/^[[:xdigit:]]+ <.*>:/ { owner = $0 } $0 ~ pattern { print owner }' \
  | sort -u > "${sse_above_baseline_owners}"
objdump -dC -- "${binary}" \
  | sed -n '/^[[:xdigit:]][[:xdigit:]]* <.*>:/h; /[[:space:]]sha\(1\|256\)/{x;p;x;}' \
  | sort -u > "${sha_owners}"
objdump -dC -- "${binary}" \
  | sed -n '/^[[:xdigit:]][[:xdigit:]]* <.*>:/h; /[[:space:]]xgetbv/{x;p;x;}' \
  | sort -u > "${xgetbv_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *bsbit_index::simd::popcnt_u64* \
      | *bsbit_align::verification::narrow::*_sse41* \
      | *bsbit_align::verification::narrow::*_sse42* \
      | *bsbit_align::verification::narrow::*_avx2* \
      | *bsbit_align::verification::narrow::*_avx512* \
      | *rans_*_sse4* \
      | *rans_*_avx2*) ;;
    *) die "POPCNT escaped a guarded backend: ${owner}" ;;
  esac
done < "${popcnt_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *bsbit_align::verification::narrow::*_avx2* \
      | *bsbit_align::verification::narrow::*_avx512* \
      | *rans_*_avx2* \
      | *rot32_simd*) ;;
    *) die "AVX escaped a guarded backend: ${owner}" ;;
  esac
done < "${avx_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *bsbit_align::verification::narrow::*_avx512*) ;;
    *) die "AVX-512 escaped its guarded backend: ${owner}" ;;
  esac
done < "${avx512_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *bsbit_align::verification::narrow::*_sse41* \
      | *bsbit_align::verification::narrow::*_sse42* \
      | *bsbit_align::verification::narrow::*_avx2* \
      | *bsbit_align::verification::narrow::*_avx512* \
      | *nibble2base_ssse3* \
      | *rans_*_sse4* \
      | *rans_*_avx2* \
      | *rot32_simd* \
      | *sha2::sha256::x86_sha::compress*) ;;
    *) die "SSE3/SSSE3/SSE4 escaped a guarded backend: ${owner}" ;;
  esac
done < "${sse_above_baseline_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *sha2::sha256::x86_sha::compress*) ;;
    *) die "SHA-NI escaped its cpufeatures-guarded dependency leaf: ${owner}" ;;
  esac
done < "${sha_owners}"
while IFS= read -r owner; do
  case "${owner}" in
    *core::core_arch::x86::xsave::_xgetbv* \
      | *htscodecs_tls_cpu_init* \
      | *__cpu_indicator_init*) ;;
    *) die "XGETBV escaped a CPUID-guarded detector: ${owner}" ;;
  esac
done < "${xgetbv_owners}"
readonly sha_instruction_count="$(
  objdump -d -- "${binary}" | rg -c '[[:space:]]sha(1|256)' || printf '0\n'
)"
readonly xgetbv_instruction_count="$(
  objdump -d -- "${binary}" | rg -c '[[:space:]]xgetbv' || printf '0\n'
)"
readonly libsais_symbol_count="$(
  objdump -t -- "${binary}" \
    | rg -c '[[:space:]]libsais(64)?\.c$' \
    || true
)"
[[ ${libsais_symbol_count} =~ ^[1-9][0-9]*$ ]] \
  || die 'binary contains no linked libsais object symbol'
printf 'rust_isa_baseline=x86-64-sse2\nruntime_dispatch=avx512f+avx512bw+avx2+popcnt,avx2+popcnt,sse4.2+popcnt,sse2,scalar\npopcnt_instruction_sites=%s\navx_instruction_sites=%s\nsse_above_baseline_instruction_sites=%s\nsha_ni_instruction_sites=%s\nxgetbv_detector_sites=%s\nwhole_elf_avx512_multiversion_sites=%s\nuncontracted_instruction_sites=%s\nlibsais_object_symbols=%s\n' \
  "${popcnt_count}" "${avx_count}" "${sse_above_baseline_count}" "${sha_instruction_count}" "${xgetbv_instruction_count}" "${avx512_count}" "${uncontracted_instruction_count}" "${libsais_symbol_count}" \
  > "${audit_dir}/instruction-contract.txt"
else
  readonly neon_instruction_count="$(
    objdump -d -- "${binary}" \
      | rg -c '(^|[[:space:],])v[0-9]+[.](8b|16b|4h|8h|2s|4s|1d|2d)([[:space:],]|$)' \
      || true
  )"
  [[ ${neon_instruction_count} =~ ^[1-9][0-9]*$ ]] \
    || die 'AArch64 binary contains no NEON instruction'
  readonly neon_backend_symbols="${audit_dir}/neon-backend-symbols.txt"
  objdump -tC -- "${binary}" \
    | rg 'bsbit_(align::verification::narrow::(neon::|.*_neon_dispatch)|index::simd::neon_popcount)' \
    > "${neon_backend_symbols}" \
    || true
  readonly neon_backend_symbol_count="$(rg -c '.' "${neon_backend_symbols}" || true)"
  [[ ${neon_backend_symbol_count} =~ ^[1-9][0-9]*$ ]] \
    || die 'AArch64 binary contains no retained bsbit NEON backend symbol'
  rg -q 'bsbit_align::verification::narrow::' "${neon_backend_symbols}" \
    || die 'AArch64 binary contains no retained alignment NEON backend symbol'
  rg -q 'bsbit_index::simd::neon_popcount' "${neon_backend_symbols}" \
    || die 'AArch64 binary contains no retained index NEON backend symbol'
  readonly sve_sme_pattern='(^|[[:space:],])(z[0-9]+[.][bhsdq]|p[0-9]+([.][bhsd])?([/][zm])?|za([0-9]+)?[.][bhsdq]|zt0)([[:space:],/]|$)|[[:space:]](smstart|smstop|rdvl|addvl|subvl|cnt[bhwd]|inc[bhwd]|dec[bhwd]|ptrue|pfalse|while[a-z]+)[[:space:]]'
  readonly sve_sme_instruction_count="$(
    objdump -d -- "${binary}" | rg -c "${sve_sme_pattern}" || printf '0\n'
  )"
  [[ ${sve_sme_instruction_count} == 0 ]] \
    || die 'AArch64 binary unexpectedly contains SVE or SME instructions'
  readonly libsais_symbol_count="$(
    objdump -t -- "${binary}" \
      | rg -c '[[:space:]]libsais(64)?\.c$' \
      || true
  )"
  [[ ${libsais_symbol_count} =~ ^[1-9][0-9]*$ ]] \
    || die 'binary contains no linked libsais object symbol'
  printf 'rust_isa_baseline=aarch64-neon\nruntime_dispatch=neon,scalar\nneon_instruction_sites=%s\nneon_backend_symbols=%s\nsve_sme_instruction_sites=%s\nlibsais_object_symbols=%s\n' \
    "${neon_instruction_count}" "${neon_backend_symbol_count}" "${sve_sme_instruction_count}" "${libsais_symbol_count}" \
    > "${audit_dir}/instruction-contract.txt"
fi
"${binary}" cpu > "${audit_dir}/runtime-cpu.txt"
if [[ ${architecture} == x86_64 ]]; then
  rg -qx 'architecture=x86_64' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report has the wrong architecture'
  rg -qx 'sse2=1' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report lost the x86-64 SSE2 baseline'
  rg -qx 'backend=(avx512|avx2|sse4[.]2|sse2|scalar)' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report selected an invalid x86-64 backend'
else
  rg -qx 'architecture=aarch64' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report has the wrong architecture'
  rg -qx 'neon=1' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report lost the AArch64 NEON baseline'
  rg -qx 'backend=(neon|scalar)' "${audit_dir}/runtime-cpu.txt" \
    || die 'built CPU report did not select NEON'
fi
sha256sum -- "${binary}" > "${audit_dir}/binary.sha256"
python3 "${repository_root}/scripts/check-release-notices.py" \
  --artifact binary --assemble "${output_root}/licenses" \
  > "${audit_dir}/license-assembly.txt"
[[ -f ${output_root}/licenses/third-party/libsais-2.10.4-LICENSE ]] \
  || die 'binary license assembly omitted libsais'
rg -q '"name"[[:space:]]*:[[:space:]]*"libsais"' \
  "${output_root}/licenses/license-manifest.json" \
  || die 'binary license manifest omitted libsais'
if [[ ${native_source_identity} == verified ]]; then
  "${script_dir}/check-native-sources.sh" \
    > "${audit_dir}/native-sources-after.txt"
  cmp -- \
    "${audit_dir}/native-sources-before.txt" \
    "${audit_dir}/native-sources-after.txt" \
    || die 'pinned native sources changed during build'
else
  printf '%s\n' \
    'unavailable: pinned native-source identity was not established before build' \
    > "${audit_dir}/native-sources-after.txt"
fi
if [[ ${git_metadata} == available ]]; then
  git -C "${repository_root}" status --short \
    > "${audit_dir}/git-status-after-build.txt"
  cmp -- \
    "${audit_dir}/git-status-before-build.txt" \
    "${audit_dir}/git-status-after-build.txt" \
    || die 'working tree changed during build'
else
  printf '%s\n' 'unavailable: Git worktree metadata could not be resolved' \
    > "${audit_dir}/git-status-after-build.txt"
fi
if (( verify_current )); then
  printf 'verified current-source binary built in %s (NOT FOR RELEASE)\n' \
    "${output_root}"
else
  printf 'audited bsbit release binary built in %s\n' "${output_root}"
fi
