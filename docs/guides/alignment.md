# Align reads

`bsbit align` maps bisulfite sequencing reads to a reference genome and writes
an input-order BAM. Read 1 alone selects single-end alignment; synchronized
read 1 and read 2 select paired-end alignment. Either layout supports
directional and non-directional libraries through the same command.

## Inputs

Alignment requires a [reference index](indexing.md) and one or two FASTQ files.
FASTQ may be plain, gzip-compressed, or BGZF-compressed. See [Input
data](../reference/input-data.md) for the complete input requirements.

## Run alignment

For paired-end data, supply both read files:

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
|---|---|---|---|
| `-i`,<br>`--index` | `PATH` | Required | Index created by `bsbit index` |
| `-1`,<br>`--read1` | `PATH` | Required | Single-end FASTQ or paired-end read 1 |
| `-2`,<br>`--read2` | `PATH` | None | Paired-end read 2 |
| `-o`,<br>`--output` | `PATH` | Required | Input-order BAM destination |
| `-t`,<br>`--threads` | `N` | `1` | Mapping workers, 1–64 |
| `--compression-threads` | `N` | `1` | BGZF output workers; 0 uses synchronous compression |
| `--compression-level` | `default\|0..9` | `1` | BGZF compression level |
| `--output-contract` | `minimal\|bismark` | `minimal` | Emit `NM/XG`, or add Bismark-compatible `MD/XM/XR` tags |
| `--mapped-only` | — | Off | Omit records without an accepted placement; retained MAPQ-0 placements remain |
| `--metrics` | — | Off | Write the layout-specific profiling TSV to standard output |

Existing regular-file destinations are atomically replaced only after the new
BAM is complete. Directories and special files are rejected.

## Configure alignment

### Search sensitivity

Default and `--sensitive` are the two public search modes for either layout.
Default mode provides the normal latency/accuracy balance. Sensitive mode
examines the wider bounded candidate frontier before final classification:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample.fastq.gz \
  -o sample.bam \
  -t 8 \
  --sensitive
```

For single-end input, sensitive mode preserves the exact default result as its
incumbent, replays seed intervals admitted by the wider 4,096-hit bound, and
completes the six-round bounded frontier. Weak, distance-three, and unmapped
incumbents retain distance-5 verification. A MAPQ-20-or-better incumbent at
edit distance two or less uses the sufficient distance-3 audit boundary. A
different-origin replacement or new rescue must be unique at MAPQ 20 or above;
a lower-confidence conflict retains the incumbent at MAPQ 0, and an
uncertified rescue remains unmapped. Pair geometry and mate rescue do not apply
to a single read.

The paired-end modes use these stable strategies:

| Mode | Stable strategy | Main behavior |
|---|---|---|
| Default | `balanced-d5-adapter-recovery-read-complete-v2` | Balanced distance-3-to-5 search plus exact adapter recovery for otherwise-unmapped directional pairs |
| `--sensitive` | `sensitive-bounded-integrated-read-complete-v1` | Additional complete-frontier, adapter, confidence, and bounded-repeat evidence |

Both modes normally emit one primary record for every input read. Ambiguous
placements can retain one deterministic representative at MAPQ 0. Use
pair-minimum MAPQ 20 when downstream analysis requires native-unique paired
evidence. bsbit MAPQ tiers rank confidence within bsbit and are not
probability-matched to another aligner. See the [performance
evidence](../performance-evidence.md) for qualified workloads and limitations.

### Directional adapter recovery

Directional single-end default and sensitive output inspect the final 30-base
domain for an exact supported prefix of the Illumina universal adapter.
Recovery requires at least 8 adapter bases and 50 retained read bases. A
tentative unique placement must remain unique at the same strand-aware origin
after 8 additional retained bases are removed.

An accepted recovery preserves the complete input sequence and qualities in
the BAM and represents the adapter tail as a strand-correct soft clip. Its MAPQ
is the lower of recovery and stability evidence, capped at 20. For a read that
was already mapped, recovery may correct the endpoint only when it preserves
the same biological origin. This is a narrow alignment fallback, not a general
adapter-trimming pipeline.

### Library type

Directional alignment searches OT and OB by default. Add `--non-directional`
to also search CTOT and CTOB and make one decision across all four supported
directions:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  -t 8 \
  --non-directional
```

For single-end input, an equal-best cross-pass result is ambiguous at MAPQ 0;
a weaker cross-pass result contributes to confidence separation and repeat
pressure. The option does not infer orientation from assay name and does not
silently approximate PBAT.

### Output contract

The default `minimal` contract writes the tags required by bsbit callers. Use
`--output-contract bismark` when a downstream consumer also requires
Bismark-style optional tags. `--output-contract` and `--mapped-only` apply to
both single-end and paired-end input.

## Paired-end resource controls

Template span, total-thread budgeting, and paired batching controls are
paired-end-only and fail explicitly for single-end input. On a dedicated host,
`--total-threads 14` selects the qualified throughput split and replaces both
explicit thread flags:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  --total-threads 14
```

A balanced stride-16 index uses 11 mapping plus 3 BGZF workers for this budget;
a fast stride-8 index uses 10 plus 4 because faster mapping exposes more output
pressure. Other budgets reserve approximately one fifth of physical cores for
output, bounded to four BGZF workers; fast indexes add one output worker from
budgets of at least 12 when possible. Keep explicit `--threads` and
`--compression-threads` when a scheduler or benchmark requires a fixed split.

See the [`bsbit align` CLI reference](../reference/cli.md#bsbit-align) for all
paired batch, queue, and template-span options.

## Profiling metrics

Normal runs keep standard output clean. With `--metrics`, single-end rows start
with `bsbit-single-alignment-metrics-v2` and paired-end rows start with
`bsbit-alignment-metrics-v2`. These counters describe search and output work;
they are not alignment-quality scores. Redirect standard output when retaining
them:

```bash
bsbit align -i GRCh38.bsbit -1 sample.fastq.gz -o sample.bam \
  --metrics > sample.align.metrics.tsv
```

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
