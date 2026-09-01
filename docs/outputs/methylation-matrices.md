# Build matrix

`bsbit combine` streams one or more sorted 18-column bsbit bedMethyl files
into sample-by-site methylation-level, count, or paired level/count matrices.
It preserves per-sample evidence and does not infer a new biological state.

## Before you run

Generate each sample with [`bsbit call meth --format bed`](methylation.md).
Inputs may be plain, gzip, or BGZF and must follow the same coordinate/context
contract. Review all combine defaults in the [CLI
reference](../reference/cli.md#bsbit-combine).

## Run combine

```bash
bsbit combine \
  --input tumor.bed.gz,normal.bed.gz,control.bed.gz \
  --sample-name tumor,normal,control \
  --output cohort.bed.gz \
  --matrix both \
  --min-count 10 \
  --min-prop 0.8 \
  --compress true \
  --threads 8
```

## Name inputs and outputs

| Choice | Rule |
|---|---|
| `--input` | Accepts comma-separated paths and may be repeated; declaration order sets matrix order |
| `--sample-name` | Supply one unique, nonempty label per input; otherwise the exact path text is used |
| Commas | Delimit values and cannot occur inside a path or label |
| `--matrix both` | Writes separate level and count files; the unsuffixed output path is not created |

The example creates `cohort.level.bed.gz` and `cohort.count.bed.gz` as one
publication operation. The [CLI reference](../reference/cli.md#bsbit-combine)
defines the complete suffix-insertion rules for every supported output name.

## Interpret filtering and output

`--min-count` filters each sample cell. A cell is valid when methylated plus
unmethylated coverage reaches the threshold. Missing sites and cells below the
threshold are `.`, never numeric zero.

`--min-prop` filters each site. With `S` samples, at least
`ceil(min_prop * S)` cells must be valid; at least one valid sample is always
required. Defaults are `--min-count 1 --min-prop 0`.

Every output begins with three `##bsbit_*` metadata lines and one `#chrom`
header. The first six columns are BED6:

| Column | Meaning |
|---:|---|
| 1–3 | Contig and 0-based half-open one-base interval |
| 4 | Modification and context, such as `m,CG,0` |
| 5 | BED score, currently `0` |
| 6 | Cytosine strand, `+` or `-` |

Remaining columns depend on the selected matrix:

| Mode | Per-sample columns | Missing value |
|---|---|---|
| `level` | `SAMPLE`, methylated / total to six decimals | `.` |
| `count` | `SAMPLE_meth_count`, `SAMPLE_total_count` | `.`, `.` |
| `both` | Separate files using the two schemas | As above |

Count output preserves methylated and total counts for downstream count-aware
models. Unmethylated count is `total_count - meth_count`. Opposite strands of
one CpG dyad remain separate rows.

## Limits and publication

Input validation covers the 18-column schema, one-base coordinates, strand,
context, counts, duplicates, and strict sorting. Conflicting context metadata
or contradictory contig order fails instead of being silently reconciled.

Memory scales mainly with samples and contigs, not total genomic sites. The
second pass opens one stream per sample, so very large cohorts may require a
higher file-descriptor limit. `--threads` controls the bounded hierarchical
merge without changing output order.

Destinations are create-only. With `both`, either both files publish or neither
does. `--compress true` writes deterministic BGZF compatible with
`tabix -p bed`.

## Next

- [Review combine parameters and suffix rules](../reference/cli.md#bsbit-combine)
- [Interpret methylation evidence](methylation.md)
- [Check supported uses and limitations](../known-differences.md)
