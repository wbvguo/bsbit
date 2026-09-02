# Rust product modules

The root Cargo workspace is split by responsibility:

| Crate | Responsibility |
|---|---|
| `bsbit-core` | Stable DNA, bisulfite-chemistry, coordinate, structural CIGAR, and reference-identity values |
| `bsbit-io` | Format-neutral file staging, identity validation, synchronization, atomic replacement, and rollback |
| `bsbit-hts` | FASTA/FASTQ/BED3/SAM/BAM and shared biological text formats, including the private audited HTSlib adapter |
| `bsbit-index` | Reference ownership, FM/rank/locate representations, `build/` construction, and `storage/` persistence formats |
| `bsbit-align` | Complete read alignment: edit distance, CIGAR/DP/Myers/SIMD verification, seeding, candidates, extension, paired placement, rescue, and MAPQ |
| `bsbit-call` | Shared fragment-evidence analysis plus independent methylation, SNP, and joint calling APIs |
| `bsbit-combine` | Bounded-memory parallel k-way merge of CGmap and/or extended bedMethyl files into count/level matrices |
| `bsbit-cli` | The single user-facing `bsbit` executable and cross-crate command orchestration |

Crate-local public API tests live in `crates/*/tests/*.rs`; checks that require
private implementation details live beside the implementation as unit tests or
in `tests/whitebox/`, and ignored local-data checks live in
`tests/qualification/`. Reusable test-only oracles shared by a crate's
integration tests live in `tests/support/`; they are not product modules.
New algorithms, feature switches, profiles, and ablations begin in an ignored
dated `workspace/worktree/` attempt; crate manifests and `src/` receive them only
after qualification and promotion under clear domain names. Shared,
small fixtures live under the workspace `tests/fixtures/`. Unsafe Rust is
confined to the exact implementation modules that own CPU intrinsics, platform
file identity, or native FFI; sibling modules and the rest of the workspace
deny it.

The production dependency graph is intentionally one-way: `io` is
format-neutral; `hts` may use `io` and `core` but never `align` or `index`;
`index` may use `core` and `io`; `align` may use `core` and `index`; `call` and
`combine` consume biological formats through `hts` and format-neutral path
validation through `io`; and `cli` performs cross-crate composition. The
policy test in `tests/tools/test_crate_boundaries.py` rejects accidental
reverse dependencies, retired experiment module names, and product crates
without a crate-level contract test.

Command modules mirror the public command names. `command/align/mod.rs` owns
shared parsing and selects single-end or paired-end layout from `--read1` and
`--read2`; parallel `single.rs` and `paired.rs` children own their FASTQ/BAM
runtimes. The mapping implementation lives in `bsbit-align`. Shared domain
identities such as `ReferenceSemanticDigest` are imported directly from
`bsbit-core`; downstream crates do not re-export them through compatibility
facades.
