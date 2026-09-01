# Testing and release checks

## Default clean-checkout gate

The source compatibility contract is Rust 1.89 or newer. CI checks the exact
1.89 lower bound and current stable independently; ordinary development may
use any toolchain in that range. The audited release artifact remains pinned
to Rust 1.94 because its standard-library license inventory is checksum-bound
to that toolchain.

Run the default gate on a supported Linux toolchain:

```sh
git submodule update --init --recursive
scripts/check-native-sources.sh
cargo fmt --all -- --check
cargo test --locked --workspace
cargo clippy --locked --workspace --all-targets --all-features \
  --no-deps -- -D warnings
cargo build --locked --release -p bsbit-cli --bin bsbit
```

This gate must pass when `agent/`, `workspace/`, `target/`, `build/`, and
`dist/` are absent.

## Feature matrix

The supported product surface is the explicit feature set in the release
build above. The live feature graph must also remain composable, so every
change runs both the release build and the workspace compatibility gate:

```sh
cargo test --locked --workspace --all-features --all-targets
cargo test --locked -p bsbit-index --no-default-features \
  --features combined-index --all-targets
cargo test --locked -p bsbit-index --no-default-features \
  --features index-construction --all-targets
cargo test --locked -p bsbit-align --all-targets
cargo test --locked -p bsbit-cli --all-targets
```

The feature surface contains product capabilities, not implementation stages:
`bsbit-index/combined-index` is the mmap reader/query closure, while
`bsbit-index/index-construction` adds the one current builder used only behind
`bsbit index`. `bsbit-align` and `bsbit-cli` have no feature switches: their
ordinary library and product surfaces are always compiled, while the CLI
enables index construction for its `bsbit index` command. The selected rank,
locate, verification,
endpoint, adapter, MAPQ, worker, and SIMD strategies are ordinary production
code behind those umbrellas. A new algorithm, mapping mode, or switch is not a
Cargo feature while it is under development: its source, runner, and evidence
stay in an ignored dated `agent/worktree/` attempt or detached worktree. The
status and ownership rules are in the [Cargo feature lifecycle](development/feature-lifecycle.md).

CI additionally links the actual `combined-index` reader/query closure and the
current alignment library while `CC` and `AR` point to nonexistent programs.
This proves alignment does not pull in index construction. The sole
`index-construction` closure still uses the pinned libsais/OpenMP toolchain.

Changes to production candidate search or classification must also run the two
mode-aware suites explicitly:

```sh
CARGO_INCREMENTAL=0 RUSTFLAGS='-A dead-code' \
  cargo test --locked -p bsbit-align
CARGO_INCREMENTAL=0 RUSTFLAGS='-A dead-code' \
  cargo test --locked -p bsbit-cli
```

The CLI test above verifies that only default and `--sensitive`
mapping modes are accepted. Development candidates must exercise the released
boundary from their dated attempt; they must not add hidden product CLI modes.

The small, standard-library-only Python policy checks run independently of
large truth data:

```sh
python3 -m unittest discover -s tests/tools -p 'test_*.py' -v
```

The crate-boundary check fixes the eight production crates, their allowed
normal dependency edges, the presence of crate-level contract tests, and the
exact supported feature inventory and expansions. Feature-bearing internal
dependencies must disable defaults so a future downstream default cannot
silently widen a product closure. The check also lexes Rust `cfg` and
`cfg_attr` predicates throughout each crate and rejects every feature
reference not declared by that crate's manifest, while ignoring comments and
string contents. Experimental features and retired module names stay outside
production source. Development-only dependencies may construct fixtures but
do not relax the production graph.

A production MAPQ policy change additionally requires the large-corpus
procedure in [performance evidence](performance-evidence.md). The maintained
small evaluator test checks tied-score admission, missing-pair recall
denominators, operating-point F1, and truth-ledger cardinality; it does not
qualify any MAPQ value.

Production search changes preserve deterministic BAM hashes for unchanged
public strategies. Sensitive release qualification additionally records at
least three same-resource
5M-pair runs and reports native unique/ambiguous/unmapped, mappability,
serialized proper pairs, all-reported accuracy, and pair-minimum MAPQ >= 10,
20, 30, and 40 subsets, with exact precision, recall, and F1 at each threshold;
no one number may stand in for the others. Both public default and sensitive
read-complete outputs are checked. The unfiltered primary-record count must
equal the input-read count. A MAPQ-only change must also compare full SAM
records with field 5 removed and obtain identical hashes before truth metrics
are accepted. Candidate comparisons and raw qualification output remain in
their ignored `agent/worktree/` attempt.

Publication tables additionally require, for each layout,
`unique + ambiguous + unmapped_or_other = total input units`, exact and
within-5-bp precision/recall at Q10/Q20/Q30/Q40, and `selected = correct +
errors` at every threshold. Primary step PR-AUC uses only Q60..Q1
score-bearing-unique points while retaining the full truth denominator. MAPQ
255 is unavailable and must produce N/A numeric-threshold and PR-AUC cells.

Before widening search or changing ranking, create a dated `agent/worktree/`
attempt and keep its strategy-oracle runner, candidate ledgers, and raw output
there. The externally visible promotion boundary is defined in the
[behavior contract](behavior-contract.md); current immutable strategy IDs and
measurements live in [performance evidence](performance-evidence.md).

Calling and matrix changes are covered by the default workspace gate. They can
also be isolated during development:

```sh
cargo test --locked -p bsbit-call
cargo test --locked -p bsbit-combine
cargo test --locked -p bsbit-hts --test bam
cargo test --locked -p bsbit-cli --test cli
```

The maintained caller tests cover authoritative FASTA projection with `MD`
ignored, required-XG plus FLAG strand handling, overlap collapse, single-sample header rules,
region-union planning, deterministic plain/BGZF publication, and VCF, CGmap,
and bedMethyl schemas. Matrix tests cover strict input validation, missing-cell and
threshold semantics, common contig-order derivation, bounded parallel merge,
and rollback-protected `both` publication. Alignment tests also cover the
complete non-directional four-strand decision and both BAM tag contracts.

## Formal extended checks

- Stable CI formats and links every target in the independent `tests/fuzz`
  workspace with its own lockfile before any long-running fuzz campaign;
- `scripts/run-rust-fuzz.sh`: Rust parser/snapshot coverage-guided fuzzing;
- `scripts/run-native-fuzz.sh`: HTSlib shim and native parser fuzzing;
- `scripts/check-htslib-shim.sh`: C ABI, sanitizer, mutation, and fault checks;
- `scripts/check-platform-publication.sh`: ext4 publication and WSL 9p
  fail-closed behavior;
- `scripts/run-release-soak.sh`: extended native and Rust mutation soak;
- `scripts/check-release-notices.py`: production dependency-license validation
  and release license-file assembly.

The release-notice lane runs both the checker's self-test and the audited
`bsbit` binary closure under exact Rust 1.94. Ordinary source compatibility remains
Rust 1.89 or newer.

Each extended runner publishes output only to an explicit ignored destination
or a private temporary directory.

`tests/tools/generate_softclip_truth.py` creates deterministic paired WGBS
truth with clean reads, adapter/low-quality 3' contamination, and deliberate
equal-best decoys. It is a development qualification input, not a runtime
dependency. The maintained checks require clean and contaminated reads to keep
the same strand-aware 5' origin, complete original reads to survive BAM
serialization, and equal-best decoys to remain ambiguous.

## Fixture policy

Small shared fixtures under `tests/fixtures/` are tracked because they are
stable, redistributable, auditable, and required by automated tests.
Crate-specific fixtures stay under `crates/*/tests/fixtures/`.

Tests are placed by the boundary they exercise:

- public API and cross-module behavior use ordinary Cargo integration tests
  under `crates/*/tests/*.rs`;
- reusable test-only oracles imported by those integration tests live under
  `crates/*/tests/support/`;
- private implementation checks live under `crates/*/tests/whitebox/` and are
  loaded only by the owning module under `cfg(test)`;
- checks that need a local index, large reference, environment variable, or
  machine-specific tool live under `crates/*/tests/qualification/` and remain
  explicitly ignored in the default gate.

Unqualified implementations are not tests and do not live under `crates/`.
Their source and runners remain in an ignored dated `agent/worktree/` attempt
until one implementation passes its predeclared promotion gate.

Large shared FASTA/FASTQ inputs, reusable human indexes, and third-party tools
belong under ignored `workspace/`. Agent-run BAMs, perf data, flamegraphs,
logs, and single-run output stay with their documented `agent/worktree/` attempt;
user-owned runs may use `workspace/benchmarks/runs/`. None is a default-test
dependency.
