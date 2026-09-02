# CLI reference

This page lists every supported `bsbit` command and option. For a task-oriented
path through the commands, start with the [Workflow](../getting-started/workflow.md).
For accepted inputs and output schemas, see [Input data](input-data.md) and
[File formats](file-formats.md).

In the **Default** column, **Required** means that the option has no default,
and **None** means that it is not set when omitted.

## Choose a command

| Command | Purpose |
|---|---|
| [`bsbit index`](#bsbit-index) | Build the alignment index |
| [`bsbit align`](#bsbit-align) | Standard single-end or paired-end alignment |
| [`bsbit call meth`](#call-meth) | Methylation calling |
| [`bsbit call snp`](#call-snp) | Bisulfite-aware diploid SNV calling |
| [`bsbit call joint`](#call-joint) | Shared methylation and SNV calling |
| [`bsbit combine`](#bsbit-combine) | Combine CGmap and/or extended bedMethyl samples into matrices |

## Help and version parameters

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-h`,<br>`--help` | — | None | Print the relevant help for any executable or subcommand and exit |
| `help` | — | None | Positional help alias for top-level `bsbit`, `bsbit call`, and `bsbit combine`; prefer `--help` in scripts |
| `-V`,<br>`--version` | — | None | Print the bsbit version and exit |

</div>

??? note "Value conventions"

    - `N` is a base-10 nonnegative integer unless a row gives a narrower range.
    - `P` is a decimal from 0 to 1 with at most nine fractional digits;
      exponent notation is not accepted. `--heterozygosity` excludes both
      endpoints.
    - `BOOL` accepts `true` or `false`; flag-only options take no value.
    - `A|B` lists the accepted literal choices.

## `bsbit index`

Build the reference index used by `bsbit align`:

```bash
bsbit index -r PATH -o PATH [-t N] \
  [--index-speed balanced|fast]
```

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-r`,<br>`--reference` | `PATH` | Required | Plain or BGZF-compressed reference genome FASTA |
| `-o`,<br>`--output` | `PATH` | Required | Generated bsbit alignment index |
| `-t`,<br>`--threads` | `N` | `1` | Number of indexing workers, 1–64 |
| `--index-speed` | `balanced\|fast` | `balanced` | Alignment-speed trade-off: `fast` uses more index storage and mapping memory |

</div>

`balanced` is recommended for most workflows. `fast` uses more storage and
alignment memory but can reduce alignment time; both modes produce the same
alignment decisions. See [Build index](../guides/indexing.md).

## `bsbit align` { #bsbit-align }

Map bisulfite sequencing reads with the index created by `bsbit index` and
write an input-order BAM. Supply only read 1 for single-end data or both read
files for paired-end data:

```bash
bsbit align \
  -x reference.bsbit \
  -1 READS_OR_R1.fastq.gz \
  [-2 R2.fastq.gz] \
  -o OUTPUT.bam \
  [OPTIONS]
```

<div class="cli-options" markdown>

| Shared option | Value | Default | Description |
|---|---|---|---|
| `-x`,<br>`--index` | `PATH` | Required | Reference index created by `bsbit index` |
| `-1`,<br>`--read1` | `PATH` | Required | Single-end FASTQ or paired-end read 1; plain, gzip, or BGZF |
| `-2`,<br>`--read2` | `PATH` | None | Paired-end read 2; plain, gzip, or BGZF |
| `-o`,<br>`--output` | `PATH` | Required | Path for the input-order BAM |
| `--sensitive` | — | Off | Search a broader set of candidate alignments |
| `--non-directional` | — | Off | Make one placement decision across all four bisulfite directions |
| `--output-contract` | `minimal\|bismark` | `minimal` | Emit `NM/XG`, or add Bismark-compatible `MD/XM/XR` tags |
| `--mapped-only` | — | Off | Omit primary records without an accepted placement; retained MAPQ-0 placements remain |
| `-t`,<br>`--threads` | `N` | `1` | Number of mapping workers, 1–64 |
| `--compression-threads` | `N` | `1` | Number of BGZF output workers; use 0 for synchronous compression |
| `--compression-level` | `default\|0..9` | `1` | HTSlib/BGZF compression setting |
| `--metrics` | — | Off | Write performance-profiling counters to standard output |

</div>

<div class="cli-options" markdown>

| Paired-only option | Value | Default | Description |
|---|---|---|---|
| `--total-threads` | `N` | None | Split one 1–64 physical-core budget between mapping and output according to index speed; conflicts with `-t` and `--compression-threads` |
| `--batch-pairs` | `N` | `16384` | Input pairs per mapping batch |
| `--alignment-queue-batches` | `N` | `2` | Bounded completed-batch queue depth |
| `--min-template-span` | `N` | `0` | Inclusive minimum template span |
| `--max-template-span` | `N` | `1000` | Inclusive maximum template span |

</div>

Paired-only options fail when read 2 is not supplied. See [Align
reads](../guides/alignment.md) for workflow guidance and [Alignment metrics
TSV](file-formats.md#alignment-metrics-tsv) for the optional profiling output.

## Calling options shared by `meth`, `snp`, and `joint`

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-i`,<br>`--input` | `PATH` | Required | Coordinate-sorted, indexed caller-compatible bsbit BAM with structured `@PG` provenance |
| `-r`,<br>`--reference` | `FASTA` | Required | Authoritative reference whose normalized content matches the BAM digest; existing FAI is used, otherwise plain FASTA is scanned in memory; BGZF requires FAI/GZI |
| `--region` | `CONTIG:START-END` | Whole dictionary | 1-based inclusive region |
| `--regions-file` | `BED` | None | Plain/gzip/BGZF BED3+ targets; conflicts with `--region` |
| `-c`,<br>`--compress` | `BOOL` | `true` | Write deterministic BGZF; `false` writes plain text |
| `-t`,<br>`--threads` | `N` | `1` | Number of regional calling workers, 1–64 |
| `--min-bq` | `N` | `20` | Minimum observed-base Phred quality, 0–93 |
| `--min-mapq` | `N` | `20` | Minimum mapping quality, 0–254 |
| `--min-depth` | `N` | `10` | Minimum qualified site depth, 1–4,294,967,295 |
| `--ignore-orphan` | — | Off | Skip paired reads without the SAM proper-pair flag; retain single-end reads |

</div>

Every module accepts one biological sample per BAM. In `joint`, these options
apply to both outputs. See [Calling BAM and
reference](input-data.md#calling-bam-and-reference) for input requirements.

## `bsbit call meth` { #call-meth }

Aggregate strand-specific methylation evidence from a bsbit BAM:

```bash
bsbit call meth \
  -i sample.prep.bam \
  -r reference.fa \
  -o sample.cgmap.gz \
  -f cgmap \
  -t 8
```

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-o`,<br>`--output` | `PATH` | Required | Output path |
| `-f`,<br>`--format` | `cgmap` or `bed` | Required | Output format: CGmap or extended bedMethyl |
| `--cg-only` | — | Off | Omit CHG and CHH sites |

</div>

See [Call methylation](../guides/methylation.md) for filtering behavior and
[CGmap](file-formats.md#cgmap-methylation-calls) or [extended
bedMethyl](file-formats.md#extended-bedmethyl-calls) for output fields.

## `bsbit call snp` { #call-snp }

Call quality-weighted, bisulfite-aware diploid SNVs:

```bash
bsbit call snp \
  -i sample.prep.bam \
  -r reference.fa \
  -o sample.vcf.gz \
  -t 8
```

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--sample-name` | `NAME` | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |
| `-o`,<br>`--output` | `PATH` | Required | VCF path |
| `--min-alt-count` | `N` | `2` | Candidate threshold and selected-ALT `LowAD` threshold, 1–4,294,967,295 |
| `--min-alt-fraction` | `P` | `0.1` | Minimum strongest-ALT candidate fraction, decimal 0–1 |
| `--min-gq` | `N` | `0` | `LowGQ` filter threshold, 0–99; 0 disables it |
| `--min-aq` | `N` | `30` | Per-ALT posterior-presence `LowAQ` threshold, 0–99 |
| `--heterozygosity` | `P` | `0.001` | Reference-divergence prior, decimal strictly between 0 and 1 |
| `--underconversion-rate` | `P` | `0.0025` | Non-conversion probability, decimal 0–1 |
| `--overconversion-rate` | `P` | `0` | Overconversion probability, decimal 0–1 |

</div>

See [Call SNVs](../guides/variant-calling.md) for BAM preparation,
the likelihood model, filters, and VCF fields.

## `bsbit call joint` { #call-joint }

Produce methylation and VCF outputs while sharing the first evidence pass:

```bash
bsbit call joint \
  -i sample.prep.bam \
  -r reference.fa \
  -m sample.cgmap.gz \
  -f cgmap \
  -v sample.vcf.gz \
  -t 8
```

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--sample-name` | `NAME` | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |
| `-m`,<br>`--meth` | `PATH` | Required | Methylation-output path |
| `-f`,<br>`--meth-format` | `cgmap` or `bed` | Required | Methylation output format: CGmap or extended bedMethyl |
| `-v`,<br>`--vcf` | `PATH` | Required | VCF output path; must differ from the methylation path |
| `--cg-only` | — | Off | Omit CHG and CHH sites from the methylation output |
| `--min-alt-count` | `N` | `2` | SNV candidate threshold and selected-ALT `LowAD` threshold, 1–4,294,967,295 |
| `--min-alt-fraction` | `P` | `0.1` | Minimum strongest-ALT candidate fraction, decimal 0–1 |
| `--min-gq` | `N` | `0` | `LowGQ` threshold, 0–99; zero disables the filter |
| `--min-aq` | `N` | `30` | Per-ALT posterior-presence `LowAQ` threshold, 0–99 |
| `--heterozygosity` | `P` | `0.001` | Reference-divergence prior, decimal strictly between 0 and 1 |
| `--underconversion-rate` | `P` | `0.0025` | Non-conversion probability, decimal 0–1 |
| `--overconversion-rate` | `P` | `0` | Overconversion probability, decimal 0–1 |

</div>

Base-quality, MAPQ, minimum-depth, region, compression, and orphan settings
apply to both outputs. `--cg-only` affects only methylation; SNV-specific
options affect only variants. The two output paths must differ. See [Call
SNVs](../guides/variant-calling.md#run-joint-calling).

## `bsbit combine`

Join named CGmap and/or extended bedMethyl samples into a BED6-plus-matrix
output:

```bash
bsbit combine \
  -i tumor.cgmap.gz,normal.cgmap.gz \
  --sample-name tumor,normal \
  -p cohort \
  -m both \
  --min-count 10 \
  --min-prop 0.8 \
  -t 8
```

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-i`,<br>`--input` | `PATH[,PATH...]` | Required, repeatable | Comma-separated sorted 8-column CGmap or 18-column bsbit extended bedMethyl inputs; formats may be mixed across samples |
| `--sample-name` | `NAME[,NAME...]` | Input paths | Optional comma-separated sample labels, one per input |
| `-p`,<br>`--prefix` | `PREFIX` | Required | Prefix for the generated matrix files |
| `-m`,<br>`--matrix` | `level\|count\|both` | `level` | Level matrix, count matrix, or both as separate files from one merge |
| `--min-count` | `N` | `1` | Per-sample minimum methylated-plus-unmethylated coverage |
| `--min-prop` | `P` | `0` | Minimum proportion of samples passing `--min-count`; at least one is always required |
| `--cg-only` | — | Off | Retain only CpG sites |
| `-c`,<br>`--compress` | `BOOL` | `true` | Write deterministic BGZF; `false` writes plain text |
| `-t`,<br>`--threads` | `N` | `1` | Number of hierarchical input-merge workers, 1–64 |

</div>

See [Build methylation matrix](../guides/methylation-matrices.md) for input
ordering, output names, filtering behavior, and matrix schemas.

## Exit behavior

The CLI uses exit code 2 for usage or unsupported-mode errors and 1 for
operational failures. Result data is written to the selected output paths;
standard output is reserved for help, version information, and optional
textual reports such as alignment metrics.
