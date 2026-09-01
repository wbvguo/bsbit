# Alignment optimization experiment and branch ledger

This page is the durable record for the local alignment optimization branches
cleaned up on 2026-09-01. It records the tested hypotheses, measured outcomes,
output-identity gates, integration point, and branch disposition so that an
experiment does not need to remain as a Git branch merely to preserve its
result.

## Canonical result and exclusions

The canonical integration branch is `codex/align-optimized-integration`. Its
qualified code commit is `def3857`; it combines the full optimization branch,
the single-end hot-path/endpoint branch, the retained metrics-v2 subset from
the throughput exploration, non-directional single-end support, and the final
single-end sensitive-accuracy audit.

`codex/single-sensitive-accuracy` and `codex/non-directional-single-end` were
explicitly excluded from the first cleanup snapshot because separate Codex
tasks were still modifying them. Once those tasks became idle and both
worktrees were clean, the non-directional tip `a237f51` was merged by
`b754ce7`, then the sensitive-accuracy tip `63d0f62` was merged by `def3857`.
Their worktrees and local branches were removed after the post-merge
qualification below. The dirty primary checkout and detached
`align-general-reference` worktree were left untouched.

## Branch disposition

| Local branch | Recorded tip | Disposition | Reason |
| --- | --- | --- | --- |
| `codex/align-full-optimization` | `5db3aa3` | delete | Fully contained in the integration branch. |
| `codex/align-performance-experiments` | `a8f4436` | delete | Fully contained in `align-full-optimization`; reports and TSV files remain in the integration branch. |
| `codex/single-end-hotpath-endpoint` | `adc3fda` | delete | Merged by integration commit `32755d7`. |
| `codex/single-end-adapter-recovery` | `dc4e07d` | delete | Fully contained in `single-end-hotpath-endpoint`. |
| `codex/paired-end-modularization` | `183d73e` | delete | Tip is the current `dev` baseline. |
| `codex/single-end-sensitive` | `6341e6f` | delete | Already an ancestor of `dev` and the integration branch. |
| `codex/alignment-locate-wavefront` | `76f599a` | delete | All source candidates regressed; only the measurements below were retained. |
| `codex/alignment-throughput` | `6ee2de0` | delete | Metrics v2 was selectively ported by `8bed5f5` and qualified by `82e4984`; the other source candidates regressed. |
| `codex/align-optimized-integration` | `def3857` | keep | Canonical qualified optimization result containing both formerly active branches. |
| `codex/single-sensitive-accuracy` | `63d0f62` | delete | Merged by `def3857`; post-merge tests and 5M validation passed. |
| `codex/non-directional-single-end` | `a237f51` | delete | Merged by `b754ce7`; post-merge tests and 5M validation passed. |

## Frozen benchmark contract

The main performance screens used the frozen five-million-fragment PE150 GRCh38
simulation; single-end runs used its five-million-read R1 file. Release builds
used `x86-64-v3`, `+popcnt`, fat LTO, and one codegen unit. Standard runs used
8 mapping plus 2 BGZF threads pinned to physical CPUs
`0,2,4,6,8,10,12,14,16,18`; thread-split experiments used at most 14 physical
cores. Candidate/control blocks were interleaved where possible, and invalid
host-noise runs remain marked in the checked-in reports rather than being
folded into medians.

Before the endpoint merge, all valid variants reproduced these frozen BAMs:

- single-end: `622a047c6f1ac4999fcc6e5afe8abc2be0d74d7bcb64acf68d9cdbf57743fc6b`;
- paired-end: `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331`.

The single-end endpoint work intentionally changed a small set of records, so
the final integration hash is recorded separately below.

## Retained experiments

| Direction | Measured result | Decision |
| --- | --- | --- |
| Bounded ungapped mismatch scan | No stable single-end wall gain; paired user CPU improved about 1--3%. | Retained as a low-risk complete-distance fast path. |
| Certified direct single-end BAM construction | Clean runs were 15.90/15.91 s versus a nearby 16.37 s bounded-only run; about 2--3% wall and 1% user-CPU improvement. | Retained with shifted-gap-tie fallback. |
| Optional stride-8 sparse SA | At 8+2: single wall -6.5% and CPU -7.5%; paired wall -10.4% and CPU -12.5%. The SA sidecar grows by 1,549,875,360 bytes (1.443 GiB) and peak RSS by about 1.44--1.48 GiB. | Retained as explicit `--index-speed fast`; stride 16 remains the compatible balanced default. |
| Native BAM error buffer | Paired BAM-write stage median 11.754 to 11.181 s (-4.9%). | Retained; buffer contents are read only on native-call failure. |
| Position-major d=3 candidate slab | Single wall -7.2% and CPU -8.3%; paired wall -3.8%, CPU -4.1%, mapping-worker CPU -6.2%. | Retained with differential kernel tests. |
| Index-aware 14-core split | Balanced index uses 11+3; fast index uses 10+4. The isolated stride-16 11+3 median was 15.22 s versus 19.73 s at 8+2 (-22.9% wall, +6.2% CPU). | Retained behind total-thread budgeting. |
| FASTQ ASCII classification LUT | Matched single-end wall -5.1% and CPU -3.0%; paired decode median -1.2%. | Retained with exhaustive byte classification tests. |
| Single-end slab/composer/endpoint work | A/B median wall 17.15 to 16.73 s (-2.45%) and CPU -2.17%, with unchanged RSS. Only 15 supported-adapter records changed in 5M reads; mapping class, MAPQ, and origin stayed unchanged while endpoint/CIGAR recovery improved. | Retained by merge commit `32755d7`. |
| Metrics v2 | Default-off BAAB: candidate mean 21.705 s / 168.395 user s versus control 21.835 s / 170.585 user s, with identical BAMs. | Retained as optional instrumentation; corrected FASTQ timer excludes queue backpressure and record-worker clocks are sampled per chunk. |

The combined full-optimization report measures the fast-index single-end build
at 13.59 s versus 16.03 s for the earlier retained stride-16 source (15.2%
faster). On paired-end with the same 14 physical cores, the conservative
combined improvement is 8--12.5%, depending on host noise and the qualified
thread/index combination.

## Rejected locate and source-shape experiments

All rows below were byte-exact. Times are candidate repeats versus the two
interleaved controls.

| Experiment | Candidate wall (s) | Control wall (s) | CPU or stage result | Decision |
| --- | ---: | ---: | --- | --- |
| Four SA rows inside paired locate | 20.85, 23.12 | 20.26, 20.75 | Mapping-worker mean 118.433 versus 112.421 s. | Reject; wider scheduling increased work/bookkeeping. |
| Four independent single-end intervals | 17.06, 17.28 | 16.86, 16.81 | Mean wall +1.99%; user CPU +2.65%. | Reject. |
| Detach prepared vectors and borrow elements | 16.87, 17.07 | 16.49, 16.65 | Mean wall +2.41%; user CPU +3.04%. | Reject. |
| Split field borrows without moving vectors | 17.16, 17.94 | 17.21, 16.43 | Mean wall +4.34%; user CPU +2.49%. | Reject. |
| Output scratch reuse | 21.31, 20.38 | 20.48, 20.36 | Candidate mean wall +2.08%. | Reject. |
| In-kernel packed d=3 loads | Stable single wall +3.3% and CPU +4--5%; paired CPU +3.3%. | N/A | Packing cost exceeded saved scalar loads. | Reject; upstream interleaving was tested separately and retained. |
| Specialized paired ungapped evaluator | 19.41 | 19.53 | Record-worker CPU 18.27 to 19.35 s. | Reject; no repeatable whole-run benefit. |
| Static d=3 substitution table | Median 13.830 | 13.815 | Wall +0.1%, CPU -0.4%. | Reject as noise-level. |

These results close the idea of exposing more rows merely by widening the
current locate wavefront. A future locate candidate must remove rank/LF work or
change the index representation, rather than schedule the same logical work in
more lanes.

## Work counters retained from the throughput exploration

The qualified directional 5M metrics-v2 run on the historical throughput
branch produced the following scientific work ledger:

| Work item | Count |
| --- | ---: |
| Directional pair passes | 5,000,000 |
| Maximal-suffix lanes | 11,455,350 |
| Maximal-suffix rank operations | 412,886,382 |
| Locate calls | 16,231,549 |
| Singleton / multi-hit locate calls | 13,638,705 / 2,592,844 |
| Located rows | 122,124,538 |
| Locate LF/rank operations | 915,815,529 |
| Emitted / locally retained candidates | 18,656,489 / 18,656,489 |
| Verified placements | 24,802,001 |
| Compatible pairs | 12,009,080 |

The observed 7.50 LF steps per row matched the stride-16 layout and was the
evidence that justified the later stride-8 experiment. A warm 10k comparison
measured 0.219 s mapping-worker time directional versus 0.600 s
non-directional (2.74x); non-directional maximal-suffix rank work was 3.59x
higher and locate/LF work 2.56x higher, while emitted candidates rose only
3.0% and verified placements only 0.07%. Most complementary-orientation work
was therefore proof work discarded before the candidate frontier.

## Final integration qualification

At `82e4984`, the integration worktree passed:

- `cargo test --workspace --all-targets`;
- `cargo clippy --workspace --all-targets --all-features -- -D warnings`;
- 26 crate-boundary tests, formatting, shell syntax, and `git diff --check`;
- full five-million-fragment single and paired output checks.

Final integration BAM SHA-256 values:

- single-end: `79c47bdd63fd4c637d054dd446dff106668700239dd0b2b44642ffb6b42584ac`;
- paired-end: `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331`.

The final paired metrics-v2 run retained the same search/locate totals and
reported 19,184,752 distinct candidate starts, 25,125,969 verified placements,
12,190,085 compatible pairs, and 7,969,678 best-pair placements.

After merging the two formerly active source branches, code commit `def3857`
again passed workspace check/tests, Clippy with warnings denied, formatting,
strict MkDocs, `git diff --check`, and all 26 crate-boundary tests. Three
stride-8 5M validation runs then produced:

| Layout | Wall (s) | BAM SHA-256 | Result |
| --- | ---: | --- | --- |
| directional single | 13.55 | `79c47bdd63fd4c637d054dd446dff106668700239dd0b2b44642ffb6b42584ac` | Exact qualified-baseline match. |
| non-directional single | 45.89 | `a0a3dd541cae1542f498482f1533bf3d95d048f939106411460b2f5a79ec08b6` | New frozen non-directional baseline. |
| directional paired | 13.95 | `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331` | Exact qualified-baseline match. |

All three BAMs passed `samtools quickcheck`; counts, flag summaries,
caller-compatible provenance, timings, RSS, input/index hashes, and complete
built-in metrics are recorded in the merged-integration evidence directory.
The two source tips were confirmed ancestors of the integration head before
their clean worktrees and local branches were removed.

## Checked-in evidence

- Exploratory report:
  `workspace/benchmark/2026-09-01--align-performance-experiments/report.md`
- Exploratory selected timings:
  `workspace/benchmark/2026-09-01--align-performance-experiments/results.tsv`
- Exploratory index measurements:
  `workspace/benchmark/2026-09-01--align-performance-experiments/index-results.tsv`
- Exploratory reproduction commands:
  `workspace/benchmark/2026-09-01--align-performance-experiments/commands.md`
- Full optimization report:
  `workspace/benchmark/2026-09-01--align-full-optimization/report.md`
- Full optimization selected timings:
  `workspace/benchmark/2026-09-01--align-full-optimization/results.tsv`
- Full optimization paired metrics:
  `workspace/benchmark/2026-09-01--align-full-optimization/paired-metrics.tsv`
- Full optimization reproduction commands:
  `workspace/benchmark/2026-09-01--align-full-optimization/commands.md`
- Merged integration report:
  `workspace/benchmark/2026-09-01--align-merged-integration/report.md`
- Merged integration timing/output summary:
  `workspace/benchmark/2026-09-01--align-merged-integration/results.tsv`
- Merged integration exact metrics:
  `workspace/benchmark/2026-09-01--align-merged-integration/single-metrics.tsv`
  and `workspace/benchmark/2026-09-01--align-merged-integration/paired-metrics.tsv`
- Merged integration reproduction commands:
  `workspace/benchmark/2026-09-01--align-merged-integration/commands.md`

The reports also name the host-local `/tmp` directories that held complete raw
runs, BAMs, binaries, indexes, and `perf.data`. Those paths are supporting
scratch evidence, not the durable record; the selected measurements and exact
commands above are committed to the repository.
