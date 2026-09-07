# Align reads

`bsbit align` maps bisulfite sequencing reads to a reference genome and writes
an input-order BAM. It supports directional and non-directional libraries with
either single-end or paired-end data.

## Inputs

Alignment requires a [reference index](indexing.md) and one or two FASTQ files.
FASTQ may be plain, gzip-compressed, or BGZF-compressed. For paired-end data,
the two files must contain matching read names in the same order and have the
same number of records. See [Input data](../input-output/input-data.md) for the
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

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-x`,<br>`--index` | Index path | Required | Reference index created by `bsbit index` |
| `-1`,<br>`--read1` | FASTQ path | Required | Single-end FASTQ or paired-end read 1; plain, gzip, or BGZF |
| `-2`,<br>`--read2` | FASTQ path | — | Paired-end read 2; plain, gzip, or BGZF |
| `-o`,<br>`--output` | BAM path | Required | Path for the input-order BAM |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of mapping workers |

</div>

## Advanced options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `--sensitive` | — | Off | Complete a broader bounded candidate and confidence frontier |
| `--non-directional` | — | Off | Align reads from a non-directional library |
| `--min-template-span` | Nonnegative integer | `0` | Minimum fragment span allowed for a concordant pair; paired-end only; `0` disables the minimum |
| `--max-template-span` | Nonnegative integer | `1000` | Maximum fragment span allowed for a concordant pair; paired-end only; also limits mate-rescue search |
| `--output-contract` | `minimal` or `bismark` | `minimal` | Add Bismark-style optional tags when required |
| `--mapped-only` | — | Off | Omit unmapped reads or read pairs; retain mapped MAPQ-0 records |
| `--max-edit-distance` | Integer from `0` to `5` | `5` | Per-read verification budget for both layouts; lower values trade sensitivity for CPU time |
| `--adapter` | `auto`, `none`, `illumina`, or sequence | `auto` | Exact 3′ adapter used as clipping evidence |
| `--adapter-min-overlap` | Integer | `8` | Minimum exact adapter overlap |
| `--adapter-max-clip` | Integer from `0` to `192` | `30` | Maximum 3′ adapter scan |
| `--soft-clip` | `auto`, `none`, or `adapter` | `auto` | Mode-qualified clipping, no clipping, or adapter-only clipping |
| `--max-soft-clip` | Integer from `0` to `192` | `30` | Maximum total clipped query bases per read |
| `--metrics` | — | Off | Write performance diagnostics to standard output |

</div>

??? note "Alignment behavior controls"

    - **Sensitivity:** The default balances speed and sensitivity. Add
      `--sensitive` to search more candidate alignments for either layout at
      the cost of additional runtime.
    - **Edit budget:** `--max-edit-distance` is primarily intended
      for controlled sensitivity/CPU trade-offs and scientific ablations. A
      lower value excludes alignments outside that whole-read edit radius and
      can change mapping rate, coordinates, and MAPQ. In paired-end runs the
      configured bound applies independently to each mate and is preserved by
      adapter and rescue remapping.
    - **Adapter and soft clipping:** `--adapter auto` and `--adapter illumina`
      both select the Illumina universal sequence. A literal A/C/G/T/N sequence
      selects a custom adapter; `none` disables adapter recognition. With
      `--soft-clip auto`, SE and default PE use exact adapter-supported repair,
      while sensitive PE may also perform candidate-local semi-global endpoint
      completion. `--soft-clip adapter` disables that generic PE completion;
      `--soft-clip none` disables every soft-clipping path. The overlap and
      clip bounds are validated before reading input and are recorded in the
      metrics TSV.
    - **Library direction:** Directional alignment is the default. Add
      `--non-directional` when required; bsbit then makes one placement
      decision across all four supported bisulfite directions. This broader
      search may take longer.
    - **Paired-end fragment span:** For paired-end input,
      `--min-template-span` and `--max-template-span` set the inclusive minimum
      and maximum outer reference span allowed for a concordant pair—the
      interval from the pair's leftmost through its rightmost mapped base. The
      defaults are `0` and `1000` bp, respectively.

      Choose these parameters to cover the library's expected fragment-size
      distribution. The bounds constrain concordant pairing, while
      `--max-template-span` also limits mate rescue. Too narrow a range can
      reject valid pairs and alter mapping status or MAPQ; too broad a range
      admits implausible pairs and increases search cost.

??? note "Performance diagnostics"

    `--metrics` replaces the human-readable progress log with structured
    performance metrics. It does not affect alignment results. See
    [Alignment metrics](../input-output/output-files.md#alignment-metrics) for
    output details.

See the [CLI reference](../reference/cli.md#bsbit-align) for option limits,
conflicts, and automatic thread allocation.

## BAM output

`bsbit align` writes BAM directly to the requested output path. The default
output is compatible with `bsbit call`.

??? note "Optional output controls"

    - `--output-contract bismark` produces Bismark-compatible output without
      changing alignment decisions.
    - `--mapped-only` omits unmapped reads or read pairs but retains mapped
      records with MAPQ 0.

## Validate the BAM

The BAM follows FASTQ input order and is not coordinate-sorted. Validate it
before continuing:

```bash
samtools quickcheck -v sample.bam
samtools flagstat sample.bam
```

## Next

- [Prepare BAM file](prepare-bam.md)
- [Alignment BAM output](../input-output/output-files.md#alignment-bam)
- [CLI reference](../reference/cli.md#bsbit-align)
- [Troubleshoot](../help/troubleshoot.md)
