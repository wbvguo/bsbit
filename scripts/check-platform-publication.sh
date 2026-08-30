#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  printf 'usage: %s EXT4_ROOT WSL_9P_ROOT\n' "$0" >&2
  exit 2
fi

readonly repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly ext4_root="$(cd -- "$1" && pwd -P)"
readonly v9fs_root="$(cd -- "$2" && pwd -P)"
readonly ext4_magic="$(stat -f -c '%t' -- "${ext4_root}")"
readonly v9fs_magic="$(stat -f -c '%t' -- "${v9fs_root}")"
readonly ext4_name_max="$(getconf NAME_MAX "${ext4_root}")"
readonly v9fs_name_max="$(getconf NAME_MAX "${v9fs_root}")"
readonly ext4_path_max="$(getconf PATH_MAX "${ext4_root}")"
readonly v9fs_path_max="$(getconf PATH_MAX "${v9fs_root}")"

if [[ "${ext4_root}" == "${v9fs_root}" ]]; then
  printf 'error: qualification roots must be distinct\n' >&2
  exit 1
fi
if [[ "${ext4_magic}" != "ef53" ]]; then
  printf 'error: expected ext4 magic ef53 at %s, observed %s\n' \
    "${ext4_root}" "${ext4_magic}" >&2
  exit 1
fi
if [[ "${v9fs_magic}" != "1021997" ]]; then
  printf 'error: expected WSL 9p magic 1021997 at %s, observed %s\n' \
    "${v9fs_root}" "${v9fs_magic}" >&2
  exit 1
fi
if [[ "${ext4_name_max}" != "255" || "${v9fs_name_max}" != "255" ]]; then
  printf 'error: exact component test requires NAME_MAX 255; ext4=%s 9p=%s\n' \
    "${ext4_name_max}" "${v9fs_name_max}" >&2
  exit 1
fi
if [[ "${ext4_path_max}" != "4096" || "${v9fs_path_max}" != "4096" ]]; then
  printf 'error: exact path test requires PATH_MAX 4096; ext4=%s 9p=%s\n' \
    "${ext4_path_max}" "${v9fs_path_max}" >&2
  exit 1
fi

ext4_temp=''
v9fs_temp=''
cleanup() {
  local path
  local status=0
  for path in "${ext4_temp}" "${v9fs_temp}"; do
    if [[ -z "${path}" ]]; then
      continue
    fi
    if [[ ! -d "${path}" ]]; then
      printf 'warning: expected private lane directory is absent: %s\n' "${path}" >&2
      status=1
      continue
    fi
    if ! rmdir -- "${path}"; then
      printf 'warning: retained nonempty platform evidence: %s\n' "${path}" >&2
      status=1
    fi
  done
  return "${status}"
}
trap cleanup EXIT

ext4_temp="$(mktemp -d "${ext4_root}/bsbit-platform-ext4.XXXXXX")"
v9fs_temp="$(mktemp -d "${v9fs_root}/.bsbit-platform-9p.XXXXXX")"
readonly ext4_temp v9fs_temp

printf 'ext4 publication lane: root=%s magic=%s name_max=%s path_max=%s temp=%s\n' \
  "${ext4_root}" "${ext4_magic}" "${ext4_name_max}" "${ext4_path_max}" "${ext4_temp}"
env TMPDIR="${ext4_temp}" TMP="${ext4_temp}" TEMP="${ext4_temp}" \
  cargo test --locked -p bsbit-io -p bsbit-hts --lib \
  --manifest-path "${repo_root}/Cargo.toml"
readonly process_boundary_test='output_component_and_read_limits_have_exact_process_boundaries'
if ! cargo test --locked --manifest-path "${repo_root}/Cargo.toml" \
  -p bsbit-cli --test cli -- --list \
  | grep -Fqx "${process_boundary_test}: test"; then
  printf 'error: platform qualification test is absent: %s\n' \
    "${process_boundary_test}" >&2
  exit 1
fi
env TMPDIR="${ext4_temp}" TMP="${ext4_temp}" TEMP="${ext4_temp}" \
  cargo test --locked --manifest-path "${repo_root}/Cargo.toml" \
  -p bsbit-cli --test cli \
  "${process_boundary_test}" -- --exact

printf 'WSL 9p fail-closed lane: root=%s magic=%s name_max=%s path_max=%s temp=%s\n' \
  "${v9fs_root}" "${v9fs_magic}" "${v9fs_name_max}" "${v9fs_path_max}" "${v9fs_temp}"
env TMPDIR="${v9fs_temp}" TMP="${v9fs_temp}" TEMP="${v9fs_temp}" \
  cargo test --locked -p bsbit-io -p bsbit-hts --lib \
  --manifest-path "${repo_root}/Cargo.toml"

printf 'Platform publication qualification OK\n'
