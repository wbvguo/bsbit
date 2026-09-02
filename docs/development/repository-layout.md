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
| `docs/` | Final MkDocs pages and compact immutable evidence referenced by them | Yes | No |
| `dev/` | Human-maintained local scripts, notebooks, experiments, source data, manuscript work, and comparison tools | No | Never |
| `workspace/` | Agent-generated attempts, notes, one-off scripts, logs, benchmark results, and generated data | No | Never |
| `target/`, `build/`, `dist/`, `artifacts/` | Regenerable output | No | Never |

`dev/` and `workspace/` are intentionally ignored in their entirety, and their
absence is normal. `dev/` is maintained by people; Agents read it by default
and modify it only when explicitly requested. `dev/data/` is the local authority
for large human-maintained source inputs. Locally, every Agent attempt uses one
informative `workspace/worktree/YYYY-MM-DD--title/` directory with a README
recording background,
design, immutable inputs/code state, results, conclusion, and follow-up. Repeated
measurements, logs, profiles, one-off agent scripts, and generated BAMs remain
with that attempt. When code isolation is required, the registered detached Git
worktree lives in that attempt's `checkout/` subdirectory and is removed with
`git worktree remove`, leaving the surrounding record intact. A top-level name
has exactly one canonical date; cross-day
activity spans belong in the README. When a successor fully replaces the same
experiment, only the latest retained result remains top-level and unique
predecessor evidence moves under its `history/` directory. Human-maintained
inputs, benchmark tooling, and third-party tools may be referenced from `dev/`
instead of copied.

Within `dev/`, `benchmarks/` contains long-lived harness and analysis code,
`data/` contains authoritative source inputs, `analysis/` contains human-owned
notebooks, `manuscript/` contains the working preprint, and `tools/` contains
local comparison-tool installations. Within `workspace/`, `worktree/` is the
dated experiment archive, `runs/` contains retained self-contained run
archives, `datasets/` contains reusable generated derivatives, and `notes/`
contains local cross-attempt synthesis. An Agent that runs a harness writes its
attempt and output under
`workspace/worktree/YYYY-MM-DD--title/`.

## Dependency direction

Product dependencies flow from the CLI into safe Rust libraries and then into
the narrow SIMD, filesystem, HTSlib, and libsais boundaries. Formal tests may
depend on product code and `tests/fixtures/`; release code must not depend on
tests, scripts, docs, `dev/`, `workspace/`, or generated output. Development
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
agent-run attempt, human-maintained local asset, or generated output.
Human-maintained scripts, notebooks, experiments, source data, manuscript work,
and tools belong under `dev/`. Raw Agent benchmark rows and complete dated
reports belong to their `workspace/worktree/` attempt; reusable internal
summaries belong under `workspace/notes/`; promoted run archives may remain
under `workspace/runs/`. Final website-facing documentation belongs in `docs/`,
while compact immutable snapshots linked by that documentation belong in
`docs/development/evidence/`; complete run archives remain local.

Within a crate, keep release implementation in `src/`, public-boundary tests in
`tests/`, checks requiring private implementation details beside the code as
unit tests or in `tests/whitebox/`, and local-data checks in
`tests/qualification/`. Test-only oracles reused by ordinary integration tests
belong in `tests/support/`; neither support code nor qualification code is a
product module. All unqualified implementation variants and their switches
belong in a dated ignored `workspace/worktree/` attempt, not under a crate. Once
qualified, promote only the selected implementation into `src/` under a stable
name and add the smallest durable test at the appropriate boundary. Once
rejected or superseded, keep it out of the live feature graph, record its
verdict and recovery commit in the tracked retired-feature registry, and keep
any full local snapshot or measurement only in `workspace/worktree/`. Git history is
the durable source archive.

A fixture moves into `tests/fixtures/` only when it is small, stable,
redistributable, auditable, and required by an automated test. Human-maintained
large FASTQ/FASTA inputs and reusable indexes remain under `dev/data/`; Agent
generated BAMs, indexes, and profiles remain under `workspace/`.
