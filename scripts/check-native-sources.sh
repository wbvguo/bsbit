#!/usr/bin/env bash
set -euo pipefail

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly repository_root="$(cd -- "${script_dir}/.." && pwd)"

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

verify_gitlink() {
  local owner=$1
  local path=$2
  local expected_commit=$3
  local mode object stage tracked_path

  read -r mode object stage tracked_path < <(
    git -C "${owner}" ls-files --stage -- "${path}" | tr '\t' ' '
  )
  [[ ${mode-} == 160000 && ${object-} == "${expected_commit}" \
    && ${stage-} == 0 && ${tracked_path-} == "${path}" ]] \
    || die "unexpected gitlink for ${path}: ${mode-} ${object-} ${stage-} ${tracked_path-}"
}

verify_submodule() {
  local label=$1
  local path=$2
  local expected_commit=$3
  local absolute_path="${repository_root}/${path}"
  local actual_commit
  local untracked

  [[ -f ${absolute_path}/.git || -d ${absolute_path}/.git ]] \
    || die "${label} is not initialized; run: git submodule update --init --recursive"
  actual_commit="$(git -C "${absolute_path}" rev-parse HEAD)"
  [[ ${actual_commit} == "${expected_commit}" ]] \
    || die "${label} commit changed: expected ${expected_commit}, found ${actual_commit}"
  # Compare content rather than relying on `git status` stat-cache output.
  # A checkout shared by Windows Git and WSL Git can report false-positive
  # worktree dirtiness even when every blob hash is unchanged.
  git -C "${absolute_path}" diff --no-ext-diff --quiet \
    --ignore-submodules=all -- . \
    || die "${label} submodule has modified tracked files"
  git -C "${absolute_path}" diff --cached --no-ext-diff --quiet \
    --ignore-submodules=all -- . \
    || die "${label} submodule has staged changes"
  untracked="$(git -C "${absolute_path}" ls-files \
    --others --exclude-standard)"
  [[ -z ${untracked} ]] \
    || die "${label} submodule has untracked files: ${untracked}"
  printf '%s source OK: commit %s\n' "${label}" "${actual_commit}"
}

command -v git >/dev/null 2>&1 || die 'missing command: git'

configured_paths="$(git -C "${repository_root}" config -f .gitmodules \
  --get-regexp '^submodule\..*\.path$' | awk '{print $2}' | LC_ALL=C sort)"
[[ ${configured_paths} == $'external/htslib\nexternal/libsais' ]] \
  || die "unexpected top-level submodule paths: ${configured_paths}"
[[ $(git -C "${repository_root}" config -f .gitmodules \
  --get submodule.external/htslib.path) == 'external/htslib' ]] \
  || die 'unexpected external/htslib submodule path'
[[ $(git -C "${repository_root}" config -f .gitmodules \
  --get submodule.external/htslib.url) == 'https://github.com/samtools/htslib.git' ]] \
  || die 'unexpected external/htslib submodule URL'
[[ $(git -C "${repository_root}" config -f .gitmodules \
  --get submodule.external/libsais.path) == 'external/libsais' ]] \
  || die 'unexpected external/libsais submodule path'
[[ $(git -C "${repository_root}" config -f .gitmodules \
  --get submodule.external/libsais.url) == 'https://github.com/IlyaGrebnov/libsais.git' ]] \
  || die 'unexpected external/libsais submodule URL'

verify_gitlink \
  "${repository_root}" \
  external/htslib \
  4b705e4fada8ee2b6b15746f725ee8ac51631803
verify_gitlink \
  "${repository_root}" \
  external/libsais \
  ce90878d784b5ff7d019300535675e4a2e22aae0

verify_submodule \
  HTSlib \
  external/htslib \
  4b705e4fada8ee2b6b15746f725ee8ac51631803
[[ $(git -C "${repository_root}/external/htslib" config -f .gitmodules \
  --get submodule.htscodecs.path) == 'htscodecs' ]] \
  || die 'unexpected htscodecs submodule path'
[[ $(git -C "${repository_root}/external/htslib" config -f .gitmodules \
  --get submodule.htscodecs.url) == 'https://github.com/samtools/htscodecs.git' ]] \
  || die 'unexpected htscodecs submodule URL'
verify_gitlink \
  "${repository_root}/external/htslib" \
  htscodecs \
  b9fc194f772e45bb0a1f44b08cbf8697a1384bae
verify_submodule \
  htscodecs \
  external/htslib/htscodecs \
  b9fc194f772e45bb0a1f44b08cbf8697a1384bae
verify_submodule \
  libsais \
  external/libsais \
  ce90878d784b5ff7d019300535675e4a2e22aae0
