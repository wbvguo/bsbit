# Build methylation matrix

`bsbit combine` merges per-sample methylation calls into a site-by-sample level
matrix, count matrix, or both. CGmap and extended bedMethyl inputs can be mixed
in the same run.

## Inputs

Use one sorted CGmap or extended bedMethyl file per sample, produced by
[`bsbit call meth`](methylation.md). Inputs may be plain, gzip-compressed, or
BGZF-compressed and must follow the same coordinate and context conventions.

## Build the matrix

```bash
bsbit combine \
  -i tumor.cgmap.gz,normal.bed.gz \
  --sample-name tumor,normal \
  -o cohort.bed.gz \
  -m both \
  --min-count 10 \
  --min-prop 0.8 \
  -c true \
  -t 8
```

## Common options

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | `PATH[,PATH...]` | Required, repeatable | Sorted CGmap or extended bedMethyl inputs; formats may be mixed |
| `--sample-name` | `NAME[,NAME...]` | Input paths | Unique sample labels in input order |
| `-o`,<br>`--output` | `PATH` | Required | Output path, or filename template when `--matrix both` is used |
| `-m`,<br>`--matrix` | `level\|count\|both` | `level` | Matrix type to write |
| `--min-count` | `N` | `1` | Minimum coverage required for each sample cell |
| `--min-prop` | `P` | `0` | Minimum proportion of samples with a valid cell at each site |
| `-c`,<br>`--compress` | `true\|false` | `false` | Whether to write BGZF-compressed output |
| `-t`,<br>`--threads` | `N` | `1` | Number of input-merge workers, from 1 to 64 |

For repeated input syntax, validation rules, and other advanced details, see
the [`bsbit combine` CLI reference](../reference/cli.md#bsbit-combine).

## Configure the matrix

**Matrix type.** `level` writes one methylated fraction per sample; `count`
writes methylated and total counts; `both` writes separate level and count
files.

**Sample names.** Input order determines matrix column order. Supply one
`--sample-name` per input, or omit the option to use the input paths as labels.

**Filtering.** `--min-count` filters individual sample cells by coverage, and
`--min-prop` filters sites by the proportion of samples with valid cells. A
missing or filtered cell is written as `.`, not numeric zero.

## Output

With `-m both`, the example writes `cohort.level.bed.gz` and
`cohort.count.bed.gz`; the unsuffixed `cohort.bed.gz` path is not created. Both
matrices use BED6 site columns followed by sample values.

See [Methylation matrices](../reference/file-formats.md#methylation-matrices)
for the level and count schemas, metadata lines, coordinates, and missing-value
representation.

## Next

- [Call methylation](methylation.md)
- [`bsbit combine` CLI reference](../reference/cli.md#bsbit-combine)
- [Methylation matrix format](../reference/file-formats.md#methylation-matrices)
