# Cargo feature lifecycle

Cargo features are supported build capabilities, not implementation-stage
switches and not an archive. Every feature must have one of the statuses
below. The release closure contains only product entry points and independently
useful library capabilities.

| Status | Meaning | Naming and placement |
|---|---|---|
| Supported capability | A product entry point or independently useful library build selected by users or release automation | Stable umbrella name; documented and covered by the release gate |
| Development candidate | A new algorithm, comparison, profile, or ablation with an owner and a current runner | Not a tracked Cargo feature; source, runner, and evidence live in a dated ignored `workspace/worktree/` attempt |
| Historical | A completed, rejected, or superseded experiment | Remove it from live Cargo features and product source; use Git history or an ignored detached worktree when recovery remains useful |
| Unknown | Inherited development code whose runner or status cannot be established | Keep it outside the tracked crate tree; resolve it in `workspace/worktree/` before promotion or retirement |

## Current release inventory

The complete tracked inventory is deliberately small:

| Crate | Supported Cargo features | Meaning |
|---|---|---|
| `bsbit-index` | `default` (empty), `combined-index`, `index-construction` | `combined-index` is the independent mmap reader/query closure; `index-construction` adds the sole current builder behind `bsbit index` |
| `bsbit-align`, `bsbit-cli`, `bsbit-core`, `bsbit-call`, `bsbit-combine`, `bsbit-hts`, `bsbit-io` | none | Their ordinary library or product surface is always compiled |

The qualified rank/locate operations, projected lanes, paired search,
verification, endpoint grouping, adapter recovery, MAPQ policy, work
scheduling, and selected SIMD kernels are implementation details. They remain
testable through their owning crate, but they are selected crate code rather
than separately composable Cargo features. `--all-features` therefore means
the union of supported product capabilities, not a synthetic mixture of
historical pipeline stages.

There is no tracked `crates/*/experiments/` incubator. Development candidates
must not be forwarded through `bsbit-cli`, compiled into a product crate, or
added to the workspace feature graph. A candidate is promoted in one change
only after its dated attempt records a deterministic runner, an acceptance
gate, the observed result, and the production tests that preserve the result.
Rejected implementations are absent from the tracked crate tree. Git history
remains their durable recovery source.

## Change checklist

Before adding or retaining a feature:

1. Create a dated `workspace/worktree/` attempt and record its owner, question, and
   predeclared promotion gate before changing a crate manifest.
2. Keep the candidate source, feature switches, runners, profiles, and outputs
   in that attempt while the result is unresolved.
3. Promote only the selected implementation into its owning crate, together
   with deterministic crate tests and an updated release closure. Do not expose
   its intermediate stages as features.
4. Verify both the explicit release build and the workspace `--all-features`
   compatibility gate; neither substitutes for the other.
5. Remove a rejected or superseded implementation from the live manifest and
   source when the conclusion is final. When a result is qualified, integrate
   the winner into the owning product capability rather than retaining an
   experiment or stage switch.
6. Put full source snapshots, runners, and measurements only in an ignored
   dated `workspace/worktree/` attempt when they have continuing local value; Git
   history remains the durable archive.
