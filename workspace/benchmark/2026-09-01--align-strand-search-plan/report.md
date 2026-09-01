# Shared strand-search plan report

## Decision

Keep the refactor. Directionality now has one owner shared by single-end and
paired-end alignment, while each layout retains its own data representation,
scoring, reduction, and output policy. This removes four-way orchestration
drift without forcing unrelated single and paired algorithms behind one large
generic interface.

`LibraryProfile` expands once, outside the candidate hot loops, to a static
conversion-pass plan:

| profile | pass plan | biological strands |
| --- | --- | --- |
| directional | original | OT, OB |
| non-directional | original + complementary | OT, OB, CTOT, CTOB |

A conversion pass owns query orientation, combined-index hit relabelling, and
the complementary paired mate permutation. Single-end still owns read-level
classification and MAPQ. Paired-end still owns template geometry, rescue,
pair scoring, pair MAPQ, and mate restoration. Cross-pass evidence is merged
only after the layout-specific pass completes.

This is deliberately not a four-lane search kernel. The already qualified
two-lane kernel and its compact workspaces remain intact, so directional work
does not pay for non-directional state and paired policy does not leak into
single-end code. Because this crate is unpublished, the follow-up cleanup
removes the obsolete `PairedLibraryProfile` alias instead of carrying a second
name for `LibraryProfile`.

The follow-up cleanup at `920aa9a` also:

- makes `PairConstraints` the sole paired profile/span input, eliminating raw
  duplicate span fields and the now-unreachable aligner span error;
- replaces paired `bool` and `0/1` mate-routing conventions with a typed
  `MateRole`, which owns lane identity, query-block orientation, and the
  canonical per-mate conversion pass;
- reuses a worker-owned primary-pass result buffer for non-directional paired
  reduction instead of allocating a clone for every mapping call;
- removes repository-unused single-end compatibility entry points; and
- derives all four layout/profile provenance modes from one exhaustive table.

## Correctness result

The parent and candidate release binaries produced byte-identical BAMs for all
four modes. Every BAM also passed `samtools quickcheck`.

| mode | SHA-256 |
| --- | --- |
| directional single | `79c47bdd63fd4c637d054dd446dff106668700239dd0b2b44642ffb6b42584ac` |
| non-directional single | `a0a3dd541cae1542f498482f1533bf3d95d048f939106411460b2f5a79ec08b6` |
| directional paired | `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331` |
| non-directional paired | `87cf57ba3ba8fbba43d851dd088493165d528085095197e8529984c43e669468` |

Internal classifications and search-work counters also match exactly between
the baseline, candidate, and cleanup binaries; the retained values are in
`single-metrics.tsv` and `paired-metrics.tsv`. A CLI integration case covers
non-directional paired dispatch and provenance, complementing the existing
directional single, non-directional single, and directional paired coverage.

## Performance result

The fixture contains five million reads or read pairs and uses the qualified
GRCh38 stride-8 index. Single-end runs used 8 mapping plus 2 BGZF workers;
paired-end runs used 10 plus 4 workers, pinned to physical cores. Both frozen
binaries used the same release flags, compression level 1, and metrics mode.

| mode | baseline wall | candidate wall | user-CPU change | interpretation |
| --- | ---: | ---: | ---: | --- |
| directional single | 14.26 s | 13.20 s | -0.55% | wall delta is dominated by a 2.42 s system-CPU change; compute is neutral |
| non-directional single | 44.14 s | 43.94 s | -0.72% | one observation; compatible with the removed result-vector allocation, but not a speedup claim |
| directional paired | 13.74 s mean | 13.785 s mean | -0.10% | two order-reversed observations; neutral |
| non-directional paired | 29.42 s | 29.66 s | +0.68% | one observation; within normal host variance |

Directional paired was repeated with the execution order reversed. Its mean
wall change is +0.33% and mean user-CPU change is -0.10%. The other modes have
one matched observation per binary, so their small deltas must not be treated
as stable estimates. Peak RSS is effectively unchanged. The defensible result
is no material throughput regression or improvement from the responsibility
refactor itself.

The responsibility refactor made single-end reuse a retained primary-pass
buffer. The cleanup applies the same ownership pattern to paired-end, removing
its per-call primary-result clone. Any benefit is below end-to-end measurement
noise and is not claimed as a throughput improvement.

## Follow-up cleanup qualification

The cleanup binary was rebuilt with the same release flags and rerun once in
each mode on the same frozen five-million-record fixture. All four BAMs are
byte-identical to the qualified `44f1d4c` candidate, all four passed
`samtools quickcheck`, and all reported classifications and search-work
counters are unchanged.

| mode | cleanup wall | cleanup user CPU | peak RSS |
| --- | ---: | ---: | ---: |
| directional single | 14.54 s | 107.04 s | 8.74 GiB |
| non-directional single | 45.18 s | 359.53 s | 8.74 GiB |
| directional paired | 12.94 s | 133.59 s | 9.35 GiB |
| non-directional paired | 28.65 s | 288.68 s | 9.35 GiB |

These are single follow-up observations taken after the original A/B campaign;
host-state and run-order effects are uncontrolled. They qualify absence of a
large regression, not a speedup. Raw values are appended to `results.tsv`.

## Quality gates

The following completed successfully:

- locked workspace check and full all-target test suite;
- `bsbit-align` and `bsbit-cli` all-target suites, including the new CLI case;
- all-target/all-feature Clippy with warnings denied;
- Rust formatting and `git diff --check`;
- 26 crate-boundary tests;
- strict MkDocs build;
- four-mode BAM byte comparison, SHA-256 capture, and `samtools quickcheck`.

Eight local frozen-fixture qualification tests remain ignored by design in the
ordinary workspace test run; the explicit five-million-record A/B campaign
above supplies the corpus-level qualification for this branch.

## Evidence identity

- baseline source: `da5ea39`
- candidate source: `44f1d4c`
- cleanup source: `920aa9a`
- baseline binary SHA-256:
  `049ff373cff4d421ba4450e7991151e9340b119c46de08ce34d0d78cb11450d6`
- candidate binary SHA-256:
  `310a2038bfb58858c7dcf568444cf86f173f43688c5ade018ccbee625dd38467`
- cleanup binary SHA-256:
  `73c3fb00bd08972be4b10a4d56dbd9b7c2126cdbea0cee30f20fe6f2d2b31c7d`
- raw binaries and BAMs: `/tmp/bsbit-align-strand-plan-20260901.FRNsk3`
- R1/R2 fixture: `/tmp/bsbit-current-benchmark-20260831/inputs`
- index: `/tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit`
- host: Intel Core i7-14700K under WSL
