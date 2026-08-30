# Pinned external source

This directory owns the pinned HTSlib and libsais Git submodules needed by a
Linux build. HTSlib recursively owns its htscodecs submodule. These are
third-party sources, not places for project experiments or generated objects.

Initialize the exact audited revisions with:

```sh
git submodule update --init --recursive
```

`scripts/check-native-sources.sh` verifies URLs, revisions, tracked content,
and untracked files. Cargo build scripts fail clearly when a submodule is not
initialized and never fetch source themselves. Project-owned wrappers live
under `crates/`; the compact production distribution-license policy lives
under `external/licenses/`.
