# Formal project scripts

This directory contains only repeatable commands needed to build or release
the current product. Scripts may create output only below an
explicit caller-supplied path, `build/`, `artifacts/`, or a temporary
directory. Those outputs are not source.

Audited release builds and formal evidence captures require committed source:
tracked, staged, or nonignored untracked changes fail before compilation
begins. `build-bsbit.sh --verify-current` runs the same compiler, CPU-report,
ISA/ELF, and license checks against an in-progress tree, but deliberately
names and marks its binary `NOT-FOR-RELEASE`; unavailable Git or native-source
identity is preserved as an audit limitation instead of being hidden.

Dated benchmark drivers, external-baseline comparison harnesses, profiling
helpers, one-off development scripts, and test-only runners do not belong
here. Test qualification, evaluation, and fuzz runners live under `tests/`.

No product binary or Cargo build may invoke this directory at runtime.

| Entry point | Purpose |
| --- | --- |
| `build-bsbit.sh` | Clean-source audited release build, or explicit `--verify-current` non-release build; runs the release-source gate before native x86-64 or AArch64 fat-LTO/optional-PGO, ISA-contract, ELF/libsais, and assembled-license checks |
| `check-release-source.sh` | Fail-fast formatting, Clippy, complete Rust test, Python policy, diff-whitespace, and duplicate-dependency gate used by formal builds |
| `check-native-sources.sh` | Verify pinned HTSlib/htscodecs and libsais submodules |
| `check-release-notices.py` | Validate the production dependency-license policy or assemble an explicit `binary` or `source` license set |
