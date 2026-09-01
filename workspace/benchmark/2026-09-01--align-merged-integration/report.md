# Merged alignment integration report

## Outcome

The merged code passed the full workspace quality gates and three real 5M
alignment runs. All BAMs passed `samtools quickcheck`, contained only primary
records, and declared the expected caller-compatible alignment mode.

The directional single-end and paired BAM SHA-256 values exactly reproduce the
previous qualified integration baselines. The non-directional single-end run
establishes a new frozen baseline for that mode.

| Layout | Input units | BAM records | Mapped records | Wall (s) | User (s) | System (s) | Peak RSS (KiB) | BAM SHA-256 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| directional single | 5,000,000 reads | 5,000,000 | 4,989,023 | 13.55 | 107.26 | 2.92 | 9,160,224 | `79c47bdd63fd4c637d054dd446dff106668700239dd0b2b44642ffb6b42584ac` |
| non-directional single | 5,000,000 reads | 5,000,000 | 4,989,028 | 45.89 | 364.69 | 4.18 | 9,162,968 | `a0a3dd541cae1542f498482f1533bf3d95d048f939106411460b2f5a79ec08b6` |
| directional paired | 5,000,000 pairs | 10,000,000 | 9,884,356 | 13.95 | 144.00 | 3.83 | 9,804,332 | `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331` |

The BAM sizes were 378,564,966 bytes, 378,579,708 bytes, and 739,799,488
bytes in the same row order. Directional single, non-directional single, and
paired provenance respectively reported
`caller-compatible-directional-single`,
`caller-compatible-nondirectional-single`, and
`caller-compatible-directional-paired`.

## Alignment metrics

Directional single classified 4,699,139 reads as unique, 289,884 as
ambiguous, and 10,977 as unmapped. It located 60,968,707 rows and verified
29,615,425 placements. Adapter endpoint recovery was attempted for 16 reads:
15 were unique, one ambiguous, none unmapped, and 135 total bases were clipped.
The direct ungapped constructor emitted 4,662,526 records; 326,497 used
traceback.

Non-directional single classified 4,698,176 reads as unique, 290,852 as
ambiguous, and 10,972 as unmapped. It located 125,087,875 rows and verified
29,625,112 placements. Adapter recovery is intentionally directional-only, so
its adapter counters are zero. It emitted 4,662,533 direct ungapped records
and 326,495 traceback records.

Paired alignment classified 4,732,194 pairs as unique, 209,984 as ambiguous,
and 57,822 as unmapped. It used the qualified stride-8 `10 mapping + 4 BGZF`
split selected by `--total-threads 14`, located 122,124,538 rows, verified
25,125,969 placements, found 12,190,085 compatible pairs, and retained
7,969,678 best-pair placements.

The single non-directional run took 3.39x the directional wall time and 4.71x
the summed mapping-worker CPU, while locating 2.05x as many rows. These are
single post-merge validation observations, not replicated performance claims.
The raw built-in rows are preserved in `single-metrics.tsv` and
`paired-metrics.tsv`.

## Quality gates

The code at `def385731ad1afb61782fe25d561e29dcaadb53a` passed:

- `cargo check --workspace --all-targets`;
- `cargo test --locked --workspace --all-targets` (eight fixture-dependent
  qualification tests remained explicitly ignored);
- `cargo clippy --locked --workspace --all-targets --all-features -- -D warnings`;
- `cargo fmt --all -- --check`;
- `git diff --check`;
- strict MkDocs build;
- all 26 crate-boundary tests.

For the three full-corpus runs, `samtools quickcheck -v` returned cleanly,
`samtools view -c` returned 5,000,000, 5,000,000, and 10,000,000 records, and
`samtools flagstat` reported zero secondary and zero supplementary records.

## Frozen inputs and environment

Release build flags were `x86-64-v3`, `+popcnt`, fat LTO, and one codegen unit.
Runs used an Intel Core i7-14700K under WSL2 Linux 6.18.33.2, Rust/Cargo 1.94.0,
and samtools/HTSlib 1.23. Single-end runs were pinned to ten physical CPUs with
8 mapping and 2 BGZF workers; paired alignment was pinned to fourteen physical
CPUs and used the automatic 10+4 stride-8 split. BAM compression level was 1.

| Artifact | SHA-256 |
| --- | --- |
| release `bsbit` | `049ff373cff4d421ba4450e7991151e9340b119c46de08ce34d0d78cb11450d6` |
| 5M R1 FASTQ | `22c4ba66773e1de2a9c6e503b1938809572f34cf11c32706c03d96d66d9c461f` |
| 5M R2 FASTQ | `0967106baf7f722c50f658840a0ec0b9b738bc47718ed4bf089544cd27adc9dd` |
| stride-8 reference catalog | `f1d2d2a876b5721f7f86c16649cce6c9432593cc610e6244b359b99d9affb53a` |
| index metadata | `290a7b84077bc916aa85bc42ec2c490678d06cf870f6c7491d0b24dcc65c6d10` |
| BWT sidecar | `4234598c76869e15ce8be48eefc5de4c7f856417a63cd1235b9fa187009960c9` |
| OCC sidecar | `3041de676132b65e599ca83adf09ce1446da1090f5523d69d4ab6c30759bd5c6` |
| SA sidecar | `a1cd1b4053f502f16c047f225829a5deb7d6df6ebc9b4060077f4ca6616bb469` |

Complete transient outputs were retained at
`/tmp/bsbit-align-final-merge-20260901.ZgPvAd` when this report was written.
