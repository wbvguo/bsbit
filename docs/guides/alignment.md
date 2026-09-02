# Align reads

`bsbit align` maps bisulfite sequencing reads to a reference genome and writes
an input-order BAM. It supports directional and non-directional libraries with
either single-end or paired-end data.

## Inputs

Alignment requires a [reference index](indexing.md) and one or two FASTQ files.
FASTQ may be plain, gzip-compressed, or BGZF-compressed. See [Input
data](../reference/input-data.md) for the complete input requirements.

## Run alignment

For paired-end data, supply both read files with `-1` and `-2`:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  -t 8
```

For single-end data, supply only read 1:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample.fastq.gz \
  -o sample.bam \
  -t 8
```

## Common options

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--index` | `PATH` | Required | Path to an alignment index created by `bsbit index` |
| `-1`,<br>`--read1` | `PATH` | Required | Plain, gzip-compressed, or BGZF-compressed FASTQ for single-end reads or paired-end read 1 |
| `-2`,<br>`--read2` | `PATH` | None | Plain, gzip-compressed, or BGZF-compressed FASTQ for paired-end read 2 |
| `-o`,<br>`--output` | `PATH` | Required | Path for the input-order BAM output |
| `-t`,<br>`--threads` | `N` | `1` | Number of mapping workers, from 1 to 64 |

For compression, batching, profiling metrics, template-span, and other
advanced options, see the [`bsbit align` CLI
reference](../reference/cli.md#bsbit-align).

## Configure alignment

### Search sensitivity

The default mode balances speed and sensitivity. Add `--sensitive` to search a
wider set of candidates when additional sensitivity is more important than
runtime:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample.fastq.gz \
  -o sample.bam \
  -t 8 \
  --sensitive
```

### Library type

Directional alignment is used by default. Add `--non-directional` for a
non-directional library; the flag works with both single-end and paired-end
data:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  -t 8 \
  --non-directional
```

### BAM output

By default, `bsbit align` writes the BAM tags required by bsbit callers. For
downstream tools that require Bismark-style optional tags, use
`--output-contract bismark`. Add `--mapped-only` to omit unmapped records.

## Validate the BAM

The BAM follows FASTQ input order and is not coordinate-sorted. Validate it
before continuing:

```bash
samtools quickcheck -v sample.bam
samtools flagstat sample.bam
```

## Next

- [Prepare BAM file](prepare-bam.md)
- [Alignment BAM fields and tags](../outputs/alignment-bam.md)
- [`bsbit align` CLI reference](../reference/cli.md#bsbit-align)
- [Troubleshoot](../help/troubleshooting.md)
