# Cross-component qualification

These test-only runners validate behavior that spans crates, native libraries,
filesystems, or CPU backends. They are not product build or runtime entry
points.

| Entry point | Purpose |
| --- | --- |
| `check-htslib-shim.sh` | Native ABI, sanitizer, mutation, and fault validation |
| `check-platform-publication.sh` | ext4 publication and WSL 9p fail-closed validation |
| `check-simd-backends.sh` | Compare every host-supported SIMD backend and require byte-identical BAM output |
| `run-release-soak.sh` | Extended release mutation and process soak |

Run these from the repository root. Each runner documents its required host
tools and arguments in its usage or initial checks.
