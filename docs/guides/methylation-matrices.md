# Build methylation matrix

`bsbit combine` combines per-sample methylation calls into a site-by-sample
matrix. It can report methylation levels, methylated and total counts, or both.

## Inputs

Supply one coordinate-sorted [CGmap](../reference/file-formats.md#cgmap-methylation-calls)
or [extended bedMethyl](../reference/file-formats.md#extended-bedmethyl-calls)
file per sample. All inputs must use the same reference genome and compatible
site coordinates. Plain, gzip-compressed, and BGZF-compressed files are
accepted.

## Build the matrix

```bash
bsbit combine \
  -i sample1.cgmap.gz,sample2.cgmap.gz \
  --sample-name sample1,sample2 \
  -p cohort \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | `PATH,...` | Required | Methylation call files, one per sample |
| `--sample-name` | `NAME,...` | Input paths | Sample names used as matrix column labels |
| `-p`,<br>`--prefix` | `PREFIX` | Required | Prefix for the generated matrix files |
| `-m`,<br>`--matrix` | `level`, `count`, or `both` | `level` | Matrix values to write |
| `-c`,<br>`--compress` | `BOOL` | `true` | Write BGZF-compressed output |
| `-t`,<br>`--threads` | `N` | `1` | Number of input-merge workers |

</div>

## Advanced parameters

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--min-count` | `N` | `1` | Minimum total coverage required to retain a sample value |
| `--min-prop` | `P` | `0` | Minimum fraction of samples that must pass `--min-count` at a site |
| `--cg-only` | — | Off | Retain only CpG sites |

</div>

??? note "Inputs and sample names"

    Separate input paths with commas. Their order determines the sample-column
    order in the matrix. CGmap and extended bedMethyl files may be mixed:

    ```bash
    -i sample1.cgmap.gz,sample2.cgmap.gz,sample3.cgmap.gz
    ```

    `--sample-name` accepts a matching comma-separated list of unique names. If
    it is omitted, the input paths are used as the column labels.

??? note "Matrix types"

    `-m level` writes one methylated fraction from 0 to 1 per sample. `-m count`
    writes methylated and total coverage for each sample. `-m both` produces
    both matrices from the same merge.

??? note "Filtering behavior"

    Add `--cg-only` to exclude CHG and CHH sites.

    `--min-count` is applied to each sample first. `--min-prop` then sets the
    fraction of samples that must pass at each site; for example, `0.8`
    requires at least 80%. Every retained site must have at least one valid
    sample, even when `--min-prop` is `0`.

See the [CLI reference](../reference/cli.md#bsbit-combine) for repeated input
syntax, accepted ranges, and complete parameter details.

## Output

The output is a coordinate-sorted BED6-plus-sample table: six genomic-site
columns followed by one methylation level per sample with `-m level`, or
methylated and total counts with `-m count`. A `.` represents a missing value
or one that does not pass `--min-count`, instead of zero.

Output names are derived from the prefix and matrix type. With `-p cohort`,
`-m level` creates `cohort.level.bed.gz`, `-m count` creates
`cohort.count.bed.gz`, and `-m both` creates both files. With `-c false`, the
files end in `.bed` instead of `.bed.gz`.

See [Methylation matrices](../reference/file-formats.md#methylation-matrices)
for the complete schemas, coordinates, and metadata fields. BGZF-compressed
output can be indexed and queried as BED:

```bash
tabix -p bed cohort.level.bed.gz
```

Index each generated file separately. The matrices can be loaded into R or
Python for sample-level quality control, clustering,
differential methylation analysis, epigenome-wide association studies (EWAS),
or methylation quantitative trait locus (mQTL) mapping.

## Next

- [Call methylation](methylation.md)
- [CLI reference: `combine`](../reference/cli.md#bsbit-combine)
- [File formats](../reference/file-formats.md#methylation-matrices)
