# Formal project scripts

This directory contains only repeatable commands needed to build, validate,
fuzz, or release the current product. Scripts may create output only below an
explicit caller-supplied path, `build/`, `artifacts/`, or a temporary
directory. Those outputs are not source.

Audited builds and formal evidence captures require committed source: tracked,
staged, or nonignored untracked changes fail before compilation begins.

Dated benchmark drivers, external-baseline comparison harnesses, profiling helpers,
and superseded milestone scripts are development history and do not belong
here. A qualification driver may remain while it is a reproducible release
gate for the current product contract.

No product binary or Cargo build may invoke this directory at runtime.

| Entry point | Purpose |
|---|---|
| `build-bsbit.sh` | Audited x86-64-v3 fat-LTO/optional-PGO umbrella build with frozen-native, ELF/libsais, and assembled-license checks |
| `check-native-sources.sh` | Verify pinned HTSlib/htscodecs and libsais submodules |
| `check-release-notices.py` | Validate the production dependency-license policy or assemble an explicit `binary` or `source` license set |
| `check-htslib-shim.sh` | Native ABI, sanitizer, mutation, and fault validation |
| `check-platform-publication.sh` | ext4 publication and WSL 9p fail-closed validation |
| `evaluate-mapq-prauc.py` | Compute tie-aware pair-level PR-AUC and Q0/Q10/Q20/Q30/Q40 operating points from a frozen truth ledger |
| `run-rust-fuzz.sh` / `summarize-rust-fuzz.py` | Rust coverage-fuzz campaign and summary |
| `run-native-fuzz.sh` / `summarize-native-fuzz.py` | Native coverage-fuzz campaign and summary |
| `run-release-soak.sh` | Extended release mutation and process soak |
