# Align reads

`bsbit align` maps bisulfite sequencing reads to a reference genome and writes
an input-order BAM. It supports directional and non-directional libraries with
either single-end or paired-end data.

## Inputs

Alignment requires a [reference index](indexing.md) and one or two FASTQ files.
FASTQ may be plain, gzip-compressed, or BGZF-compressed. For paired-end data,
the two files must contain matching read names in the same order and have the
same number of records. See [Input data](../reference/input-data.md) for the
complete input requirements.

## Run alignment

For paired-end data, supply both read files:

```bash
bsbit align \
  -x GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  -t 8
```

For single-end data, supply only read 1:

```bash
bsbit align \
  -x GRCh38.bsbit \
  -1 sample.fastq.gz \
  -o sample.bam \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-x`,<br>`--index` | `PATH` | Required | Reference index created by `bsbit index` |
| `-1`,<br>`--read1` | `PATH` | Required | Single-end FASTQ or paired-end read 1; plain, gzip, or BGZF |
| `-2`,<br>`--read2` | `PATH` | None | Paired-end read 2; plain, gzip, or BGZF |
| `-o`,<br>`--output` | `PATH` | Required | Path for the input-order BAM |
| `-t`,<br>`--threads` | `N` | `1` | Number of mapping workers, 1–64 |

</div>

## Advanced parameters

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--sensitive` | — | Off | Search a broader set of candidate alignments |
| `--non-directional` | — | Off | Align reads from a non-directional library |
| `--output-contract` | `minimal` or `bismark` | `minimal` | Add Bismark-style optional tags when required |
| `--mapped-only` | — | Off | Omit reads or read pairs without an accepted placement |
| `--compression-threads` | `N` | `1` | Number of BGZF output workers; `0` uses synchronous compression |
| `--compression-level` | `default` or `0`–`9` | `1` | BGZF compression level |
| `--metrics` | — | Off | Write performance diagnostics to standard output |

</div>

Paired-end data also accepts the following parameters:

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--total-threads` | `N` | None | Split a 1–64 core budget between mapping and output workers |
| `--batch-pairs` | `N` | `16384` | Number of read pairs per mapping batch |
| `--alignment-queue-batches` | `N` | `2` | Number of completed batches held for output |
| `--min-template-span` | `N` | `0` | Minimum accepted template span, inclusive |
| `--max-template-span` | `N` | `1000` | Maximum accepted template span, inclusive |

</div>

??? note "Sensitive alignment"

    The default mode balances speed and alignment sensitivity. Add
    `--sensitive` to search a broader set of candidate alignments. It works
    with single-end and paired-end data but may take longer:

    ```bash
    bsbit align \
      -x GRCh38.bsbit \
      -1 sample.fastq.gz \
      -o sample.bam \
      -t 8 \
      --sensitive
    ```

??? note "Non-directional libraries"

    Directional alignment is used by default. Add `--non-directional` for a
    non-directional library; bsbit then makes one placement decision across all
    four supported bisulfite directions:

    ```bash
    bsbit align \
      -x GRCh38.bsbit \
      -1 sample_R1.fastq.gz \
      -2 sample_R2.fastq.gz \
      -o sample.bam \
      -t 8 \
      --non-directional
    ```

??? note "Paired-end template span"

    Template span is the number of reference bases covered from the outer
    start of one mate to the outer end of the other. The accepted range is
    0–1000 bp by default, inclusive. Change `--min-template-span` or
    `--max-template-span` only when the expected fragment sizes require
    different bounds.

??? note "BAM output options"

    Use `--output-contract bismark` only when a downstream tool requires
    Bismark-style optional tags. `--mapped-only` omits reads or read pairs
    without an accepted placement; accepted MAPQ-0 placements remain in the
    BAM.

??? note "Thread and batching controls"

    `-t` is sufficient for most runs. `--total-threads` automatically divides
    a paired-end core budget between mapping and BAM output and cannot be used
    with `-t` or `--compression-threads`. Batch and queue settings normally do
    not need adjustment.

??? note "Performance metrics"

    `--metrics` writes runtime and workload diagnostics without changing the
    BAM. These metrics describe performance, not alignment quality. Redirect
    standard output to save them:

    ```bash
    bsbit align \
      -x GRCh38.bsbit \
      -1 sample.fastq.gz \
      -o sample.bam \
      --metrics > sample.align.metrics.tsv
    ```

See the [CLI reference](../reference/cli.md#bsbit-align) for parameter limits,
conflicts, and automatic thread allocation.

## BAM output

By default, the BAM includes the tags required by `bsbit call`. It is published
only after alignment completes successfully.

## Validate the BAM

The BAM follows FASTQ input order and is not coordinate-sorted. Validate it
before continuing:

```bash
samtools quickcheck -v sample.bam
samtools flagstat sample.bam
```

## Next

- [Prepare BAM file](prepare-bam.md)
- [Alignment BAM output](../outputs/index.md#alignment-bam)
- [CLI reference](../reference/cli.md#bsbit-align)
- [Troubleshoot](../help/troubleshooting.md)
