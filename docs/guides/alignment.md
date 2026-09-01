# Align reads

`bsbit align` is the entry point for both read layouts. Supplying read 1 alone
selects directional single-end alignment; adding synchronized read 2 selects
paired-end alignment without changing the command.

Both layouts open the opaque index created by `bsbit index`; alignment never
builds or modifies index data.

## Before you run

- Complete the [installation](../getting-started/installation.md) and make the
  built binary available on `PATH`.
- [Build one complete index](indexing.md) from the reference FASTA.
- For paired input, confirm that R1 and R2 end together and use synchronized
  read names.
- Choose template-span bounds appropriate to a paired library rather than
  relying on a learned insert-size distribution.

The audited fat-LTO build used for frozen benchmark reproduction is documented
with the [performance evidence](../performance-evidence.md#reproduction-protocol).
It is not required for an ordinary source build or analysis run.

## Align single-end reads

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample.fastq.gz \
  --output-bam sample.bam \
  --sensitive \
  --threads 8
```

Omit `--sensitive` for the low-latency default. With it, single-end alignment
replays seed intervals admitted only by the wider 4,096-hit bound, completes
the six-round bounded frontier, and then performs distance-5 verification,
origin classification, and MAPQ. Pair geometry, mate rescue, and paired
adapter recovery do not apply to a single read.

!!! note "Single-end confidence"
    Unique origins receive numeric Q10/Q15/Q20/Q30/Q40 from evidence retained
    by the selecting search; tied origins remain MAPQ 0. Output declares
    `caller-compatible-directional-single`. The published truth qualification
    is scoped to the documented 5M-R1 simulated corpus.

Single-end output preserves input order. Coordinate-based analysis still
requires coordinate sorting. Single-end input currently supports the
directional library model; non-directional single-end and PBAT are unsupported.

Shared options including `--sensitive`, `--threads`, `--bam-threads`, and
`--bam-compression-level` apply to the read-1-only layout. Paired-only controls
including `--non-directional`, template span, `--mapped-only`, output-contract
selection, paired batching, and `--metrics` fail explicitly on single-end
input instead of being ignored.

## Align paired reads

Use synchronized R1 and R2 for caller-compatible MAPQ and BAM provenance:

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --threads 8 \
  --bam-threads 2 \
  --output-contract minimal \
  --min-template-span 0 \
  --max-template-span 1000
```

The default run keeps stdout clean. Add `--metrics` and redirect stdout when
the two-row profiling TSV is useful:

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --threads 8 \
  --metrics \
  > sample.summary.tsv
```

The BAM is written in input order. Inspect it immediately, then follow the
[BAM-preparation guide](prepare-bam.md) before coordinate indexing or calling.

## Choose an alignment mode

Both layouts expose default and sensitive search. For single-end input,
default keeps the existing early-resolution d3/d5 path and `--sensitive`
completes the wider bounded candidate frontier before classification. The
paired-end path additionally records one of these stable strategies:

| Mode | Stable strategy | Main behavior |
|---|---|---|
| Default (omit a mode flag) | `balanced-d5-adapter-recovery-read-complete-v2` | Balanced distance-3-to-5 search plus exact Illumina-adapter recovery for otherwise-unmapped directional pairs |
| `--sensitive` | `sensitive-bounded-integrated-read-complete-v1` | Additional complete-frontier, adapter, confidence, and bounded-repeat evidence within the published wall-time envelope |

Use default mode for the normal latency/accuracy balance. Use `--sensitive`
when the documented extra recall is worth the additional runtime. Current
measurements and their workload boundary live on the
[performance page](../performance-evidence.md).

Both modes write one primary record per input read by default. An ambiguous
pair may retain one deterministic mapped representative, normally at MAPQ 0;
a narrowly qualified sensitive subset may receive MAPQ 10. Use pair-minimum
MAPQ 20 when a downstream workflow requires pairs classified as unique by this
executable. bsbit MAPQ tiers rank confidence within bsbit and are not
probability-matched to another aligner.

## Select library orientation

Directional paired libraries are the default. Add `--non-directional` to make
one placement decision across all four supported directional classes. The
option does not infer orientation from assay name and does not change the
output-tag contract.

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.nondirectional.bam \
  --non-directional \
  --threads 8
```

PBAT is unsupported rather than silently approximated. Check the actual
library protocol before choosing this option.

## Select the output contract

The default `minimal` contract emits `NM` and `XG` and is sufficient for the
in-tree callers. Use `--output-contract bismark` only for a consumer that
requires Bismark-compatible `MD`, `XM`, and `XR` tags in addition to `NM` and
`XG`. The compatibility contract changes tags, not coordinates, ambiguity,
MAPQ, or classification.

See [Alignment BAM](../outputs/alignment-bam.md) for header provenance, record
completeness, tag definitions, and inspection commands.

## Next

- [Prepare the BAM for calling](prepare-bam.md)
- [Review every alignment option](../reference/cli.md#bsbit-align)
- [Review supported workflows](../getting-started/workflow.md#supported-workflows)
- [Troubleshoot alignment](../help/troubleshooting.md)
