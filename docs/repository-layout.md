# Repository layout

This repository separates files by lifecycle and dependency role. A clean
checkout must build and test without local research material, large datasets,
previous build output, or benchmark results.

| Path | Responsibility | Tracked in Git | May product code depend on it? |
|---|---|---:|---:|
| `crates/*/src/` | Rust product modules | Yes | Yes |
| `crates/bsbit-hts/htslib-shim/` | Project-owned narrow C ABI over the pinned HTSlib source | Yes | Through `bsbit-hts` only |
| `crates/*/tests/` | Crate-owned public API, white-box, qualification tests, and small fixtures | Yes | Tests only |
| `external/` | Pinned third-party Git submodules plus the static production distribution-license policy under `external/licenses/` | Gitlinks and policy files | Yes, through audited wrappers only |
| `tests/` | Cross-crate fixtures, formal test support, and the independent Rust fuzz workspace | Yes | Tests only |
| `scripts/` | Reproducible build, release, and formal validation entry points | Yes | No |
| `docs/` | Final pages published by the MkDocs website | Yes | No |
| `agent/` | Agent-owned independent worktrees, new algorithm/feature incubation, and complete historical attempt records under `agent/worktree/` | No | Never |
| `workspace/` | Shared local datasets/tools, reusable or user scripts, local result summaries, manuscript work, and user-owned runs | No | Never |
| `target/`, `build/`, `dist/`, `artifacts/` | Regenerable output | No | Never |

`agent/` and `workspace/` are intentionally ignored in their entirety. Their
absence is normal. Locally, every agent attempt uses one informative
`agent/worktree/YYYY-MM-DD--topic/` directory with a README recording background,
design, immutable inputs/code state, results, conclusion, and follow-up. Repeated
measurements, logs, profiles, one-off agent scripts, and generated BAMs remain
with that attempt. When code isolation is required, the registered detached Git
worktree lives in that attempt's `checkout/` subdirectory and is removed with
`git worktree remove`, leaving the surrounding record intact. A top-level name
has exactly one canonical date; cross-day
activity spans belong in the README. When a successor fully replaces the same
experiment, only the latest retained result remains top-level and unique
predecessor evidence moves under its `history/` directory. Shared inputs and
third-party tools may be referenced from `workspace/` instead of copied.

`workspace/` is not a second experiment archive. It is the stable local supply
surface for `datasets/`, `tools/`, reusable/user-runnable `code/` and `scripts/`,
cross-attempt summaries and compact evidence under `docs/`, manuscript work,
and runs explicitly initiated by the user. An agent that runs the same harness
writes its attempt and outputs under `agent/worktree/`.

## Dependency direction

Product dependencies flow from the CLI into safe Rust libraries and then into
the narrow SIMD, filesystem, HTSlib, and libsais boundaries. Formal tests may
depend on product code and `tests/fixtures/`; release code must not depend on
tests, scripts, docs, `agent/`, `workspace/`, or generated output. Development
candidates exercise tracked code through stable public/test boundaries or a
detached worktree; they do not become dependencies of the live workspace.

Native-source checks verify the exact submodule URLs, commits, tracked bytes,
and absence of untracked files. Build scripts require initialized submodules;
they must not download dependencies or consult a developer checkout.
The release-license checker combines `external/licenses/`, the pinned native
license files, and the audited Rust toolchain copyright inventory. Assembly
requires an explicit `binary` or `source` scope and writes only
the applicable root and third-party texts into ignored `dist/`. The complete
CLI includes libsais inside `bsbit index`; the canonical `bsbit align` runtime
only reads the resulting opaque index. The binary scope also
include the Rust standard-library inventory. Product code does not read this
policy at runtime.

## Placement checklist

Before adding a file, decide whether it is required by the shipped product, a
repeatable formal test, current documentation, fixed external source, an
agent-run attempt, shared local work data/tooling, or generated output. Raw
agent benchmark rows and complete dated reports belong to their
`agent/worktree/` attempt; reusable internal summaries belong under
`workspace/docs/`; user-owned runs may remain under
`workspace/benchmarks/runs/`. Only final website-facing documentation and the
compact immutable TSV snapshots linked by that documentation belong in
`docs/`; complete run archives remain local.

Within a crate, keep release implementation in `src/`, public-boundary tests in
`tests/`, checks requiring private implementation details beside the code as
unit tests or in `tests/whitebox/`, and local-data checks in
`tests/qualification/`. Test-only oracles reused by ordinary integration tests
belong in `tests/support/`; neither support code nor qualification code is a
product module. All unqualified implementation variants and their switches
belong in a dated ignored `agent/worktree/` attempt, not under a crate. Once
qualified, promote only the selected implementation into `src/` under a stable
name and add the smallest durable test at the appropriate boundary. Once
rejected or superseded, keep it out of the live feature graph, record its
verdict and recovery commit in the tracked retired-feature registry, and keep
any full local snapshot or measurement only in `agent/worktree/`. Git history is
the durable source archive.

A fixture moves into `tests/fixtures/` only when it is small, stable,
redistributable, auditable, and required by an automated test. Large FASTQ,
FASTA, BAM, index images, and profiles remain local work data.
