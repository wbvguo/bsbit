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
  -i sample1.CGmap,sample2.CGmap \
  --sample-name sample1,sample2 \
  -o cohort.bed \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | Comma-separated file paths | Required | Methylation call files, one per sample |
| `-o`,<br>`--output` | Output path | Required | Destination; with `both`, a template for the level and count paths |
| `-c`,<br>`--compress` | `true` or `false` | `false` | Write BGZF-compressed output when `true` |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of input-merge workers |
| `--matrix` | `level`, `count`, or `both` | `level` | Matrix values to write |
| `--sample-name` | Comma-separated sample names | Input paths | Sample names used as matrix column labels |

</div>

## Advanced options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `--min-count` | Positive integer | `1` | Minimum total coverage required to retain a sample value |
| `--min-prop` | Decimal from `0` to `1` | `0` | Minimum fraction of samples that must pass `--min-count` at a site |
| `--cg-only` | — | Off | Retain only CpG sites |

</div>

??? note "Inputs and sample names"

    Separate input paths with commas. Their order determines the sample-column
    order in the matrix. CGmap and extended bedMethyl files may be mixed:

    ```bash
    -i sample1.CGmap.gz,sample2.CGmap.gz,sample3.CGmap.gz
    ```

    `--sample-name` accepts a matching comma-separated list of unique names. If
    it is omitted, the exact input path text is used as the column label.
    Commas delimit values and therefore cannot occur inside a path or label.

??? note "Matrix types"

    `--matrix level` writes one methylated fraction from 0 to 1 per sample.
    `--matrix count` writes methylated and total coverage for each sample.
    `--matrix both` produces both matrices from the same merge.

??? note "Filtering behavior"

    Add `--cg-only` to exclude CHG and CHH sites.

    `--min-count` is applied to each sample first. `--min-prop` then sets the
    fraction of samples that must pass at each site; for example, `0.8`
    requires at least 80%. Every retained site must have at least one valid
    sample, even when `--min-prop` is `0`.

See the [CLI reference](../reference/cli.md#bsbit-combine) for accepted input,
and complete option details.

## Output

The output is a coordinate-sorted BED6-plus-sample table: six genomic-site
columns followed by one methylation level per sample with `--matrix level`, or
methylated and total counts with `--matrix count`. A `.` represents a missing value
or one that does not pass `--min-count`, instead of zero.

For one matrix type, `--output` is the exact destination. With
`--output cohort.bed --matrix both`, bsbit creates
`cohort.level.bed` and `cohort.count.bed`;

See [Methylation matrices](../reference/file-formats.md#methylation-matrices)
for the complete schemas, coordinates, and metadata fields. With `-c true`, use
a `.gz` filename; BGZF-compressed output can then be indexed and queried as BED:

```bash
tabix -p bed cohort.level.bed.gz
```

Index each generated file separately. The matrices can be loaded into R or
Python for sample-level quality control, clustering,
differential methylation analysis, epigenome-wide association studies (EWAS),
or methylation quantitative trait locus (mQTL) mapping.

## Next

- [Call methylation](methylation-calling.md)
- [CLI reference: `combine`](../reference/cli.md#bsbit-combine)
- [File formats](../reference/file-formats.md#methylation-matrices)
