# CLI reference

This page lists every supported `bsbit` command and option. For a task-oriented
path through the commands, start with the [Workflow](../getting-started/workflow.md).
For accepted inputs and output schemas, see [Input data](../input-output/input-data.md) and
[File formats](file-formats.md).

## Conventions

The **Argument** column describes the value accepted after an option; **—**
marks a switch that takes no argument.

The **Default** column uses **Required** for mandatory options, **—** when no
optional input is selected, and **Off** for flags that are off unless supplied.

- Integer and decimal bounds are inclusive unless stated otherwise.
- Decimal probabilities accept at most nine fractional digits; exponent
  notation is not accepted.
- Literal choices are case-sensitive and shown in code font.

`index`, `align`, and `call` open their output paths before main input
processing. Existing regular files are overwritten and a failed command may leave
an empty or partial result. See [Filesystem and output](../help/troubleshoot.md#output-cannot-be-opened-for-writing).

## Choose a command

| Command | Purpose |
| --- | --- |
| [`bsbit index`](#bsbit-index) | Build the alignment index |
| [`bsbit align`](#bsbit-align) | Standard single-end or paired-end alignment |
| [`bsbit cpu`](#bsbit-cpu) | Report CPU features and the selected SIMD backend |
| [`bsbit call meth`](#call-meth) | Methylation calling |
| [`bsbit call snp`](#call-snp) | Bisulfite-aware diploid SNV calling |
| [`bsbit call joint`](#call-joint) | Shared methylation and SNV calling |
| [`bsbit combine`](#bsbit-combine) | Combine CGmap and/or extended bedMethyl samples into matrices |

## Help and version

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-h`,<br>`--help` | — | Off | Print the relevant help for any executable or subcommand and exit |
| `-v`,<br>`--version` | — | Off | Print the bsbit version and exit |

</div>

`help` is also accepted as a positional alias for `--help` by top-level
`bsbit`, `bsbit call`, and `bsbit combine`; prefer `--help` in scripts.

## `bsbit cpu` { #bsbit-cpu }

Report the process-visible CPU features and the backend that indexing and
alignment would select:

```bash
bsbit cpu [--simd-backend auto|scalar|sse2|sse4.2|avx2|avx512|neon]
```

The stable key/value output reports the detected architecture and CPU features,
the selected backend, and its complete instruction-set requirements. Features
that do not apply to the current architecture are reported as zero; SIMD
features and POPCNT are detected independently.

Automatic backend selection follows this priority:

- x86-64: AVX-512 → AVX2 + POPCNT → SSE4.2 + POPCNT → SSE2 → scalar
- AArch64: NEON → scalar
- Other architectures: scalar

The scalar backend can always be forced. A forced SIMD backend fails safely if
it targets the wrong architecture or requires unavailable CPU features.

??? example "Inspect CPU instruction sets on Linux"
    Run this in the same Linux, WSL2, VM, container, or batch allocation that
    will run bsbit:

    ```bash
    lscpu | grep --color=auto -E \
      'sse2|sse4_1|sse4_2|popcnt|avx2|avx512f|avx512bw|asimd|neon|$'
    ```

    The report shows the CPU architecture, logical CPU count, model, topology,
    caches, virtualization, and OS-visible instruction sets. Highlighted flags
    are those used by bsbit backends; on AArch64, NEON usually appears as
    `asimd`.

## `bsbit index`

Build the reference index used by `bsbit align`:

```bash
bsbit index -r PATH -o PATH [-t N] \
  [--memory-mib N] [--index-speed fast|compact] \
  [--simd-backend auto|scalar|sse2|sse4.2|avx2|avx512|neon]
```

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-r`,<br>`--reference` | FASTA path | Required | Plain, gzip, or BGZF reference genome FASTA |
| `-o`,<br>`--output` | Index path | Required | Path for the bsbit reference index |
| `-t`,<br>`--threads` | Positive integer | `1` | Indexing-worker count in the signed 32-bit native domain |
| `--index-speed` | `fast` or `compact` | `fast` | `fast` favors lookup speed; `compact` reduces storage and memory use |
| `--memory-mib` | Positive integer | `9300` | Search-index construction budget in MiB |
| `--simd-backend` | `auto`, `scalar`, `sse2`, `sse4.2`, `avx2`, `avx512`, or `neon` | `auto` | Detect once, or force a validated backend for qualification |

</div>

For most production workflows, use the default `fast` index. Choose `compact`
when reducing index size and peak alignment memory is more important than
alignment speed:

```bash
bsbit index \
  -r GRCh38.fa \
  -o GRCh38.compact.bsbit \
  -t 8 \
  --index-speed compact
```

On GRCh38, `compact` reduces the index size by about 1.44 GiB and lowers peak
alignment memory without changing alignment results. In benchmarks of five
million single-end reads or read pairs, `fast` completed alignment 6.5% faster
for single-end data and 10.4% faster for paired-end data than `compact`.

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

| Shared option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-x`,<br>`--index` | Index path | Required | Reference index created by `bsbit index` |
| `-1`,<br>`--read1` | FASTQ path | Required | Single-end FASTQ or paired-end read 1; plain, gzip, or BGZF |
| `-2`,<br>`--read2` | FASTQ path | — | Paired-end read 2; plain, gzip, or BGZF |
| `-o`,<br>`--output` | BAM path | Required | Path for the input-order BAM |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of mapping workers |
| `--sensitive` | — | Off | Complete the broader bounded candidate and confidence frontier for either layout |
| `--non-directional` | — | Off | Make one placement decision across all four bisulfite directions |
| `--output-contract` | `minimal` or `bismark` | `minimal` | Emit `NM/XG`, or add Bismark-compatible `MD/XM/XR` tags |
| `--mapped-only` | — | Off | Omit primary records without an accepted placement; retained MAPQ-0 placements remain |
| `--max-edit-distance` | Integer from `0` to `5` | `5` | Maximum per-read edit distance for either layout; PE applies the same bound to each mate and every recovery phase |
| `--adapter` | `auto`, `none`, `illumina`, or an A/C/G/T/N sequence | `auto` | Select exact 3′ adapter evidence; `auto` currently selects the Illumina universal adapter |
| `--adapter-min-overlap` | Integer from `1` through adapter length | `8` | Minimum exact adapter overlap required to admit an adapter boundary |
| `--adapter-max-clip` | Integer from `0` to `192` | `30` | Maximum 3′ suffix inspected for adapter evidence |
| `--soft-clip` | `auto`, `none`, or `adapter` | `auto` | Use mode-qualified clipping, disable clipping, or allow only exact adapter-supported clipping |
| `--max-soft-clip` | Integer from `0` to `192` | `30` | Maximum total query bases soft-clipped per read; `0` disables clipping |
| `--simd-backend` | `auto`, `scalar`, `sse2`, `sse4.2`, `avx2`, `avx512`, or `neon` | `auto` | Detect once, or force a validated architecture backend |
| `--total-threads` | Positive integer | — | Split one core budget between mapping and output; conflicts with both explicit thread flags |
| `--batch-size` | Positive integer | `1000` single;<br>`16384` paired | Reads or read pairs per mapping batch |
| `--queue-batches` | Positive integer | `2` | Bounded queue depth between pipeline stages |
| `--compression-level` | `default` or integer from `0` to `9` | `1` | HTSlib/BGZF compression setting |
| `--metrics` | — | Off | Suppress human progress and write a layout-specific profiling TSV to standard output |
| `--compression-threads` | Nonnegative integer | `1` | Number of BGZF output workers; use 0 for synchronous compression |

</div>

??? note "Paired-end template span"
    Template span is the number of reference bases from the pair's leftmost
    aligned base to its rightmost aligned base, including both mates and the
    interval between them. The bounds reject implausibly distant mate
    combinations and limit mate-rescue searches.

    <div class="cli-options" markdown>

    | Option | Argument | Default | Description |
    | --- | --- | --- | --- |
    | `--min-template-span` | Nonnegative integer | `0` | Minimum fragment span allowed for a concordant pair; paired-end only; `0` disables the minimum |
    | `--max-template-span` | Nonnegative integer | `1000` | Maximum fragment span allowed for a concordant pair; paired-end only; also limits mate-rescue search |

    </div>

    Both bounds are inclusive and accepted only with paired-end input. Most
    users should keep the defaults; change them only when the library's
    expected fragment sizes fall outside this range.

FASTQ sequence and quality lengths must both be 3–192 bytes. Metrics start with
`bsbit-alignment-metrics-paired-end-v1` for paired-end or
`bsbit-alignment-metrics-single-end-v1` for single-end input. See [Align
reads](../guides/alignment.md) for workflow guidance and [Alignment metrics
TSV](file-formats.md#alignment-metrics-tsv) for the optional profiling output.

## Calling options shared by `meth`, `snp`, and `joint`

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | BAM path | Required | Coordinate-sorted, BAI/CSI-indexed bsbit BAM with per-contig `@SQ M5` checksums and mapped-record `XG` tags |
| `-r`,<br>`--reference` | FASTA path | Required | Matching uncompressed or BGZF FASTA; `.fai` is recommended for uncompressed input, while BGZF requires `.fai` and `.gzi` |
| `-c`,<br>`--compress` | `true` or `false` | `false` | Write deterministic BGZF when `true`; otherwise write plain text |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of regional calling workers |
| `--region` | `CONTIG:START-END` | Whole dictionary | One 1-based inclusive region |
| `--regions-bed` | BED path | — | Plain, gzip, or BGZF BED3+ file for multiple regions |
| `--min-bq` | Integer from `0` to `93` | `20` | Minimum observed-base Phred quality |
| `--min-mapq` | Integer from `0` to `254` | `20` | Minimum mapping quality |
| `--min-depth` | Positive integer | `10` | Minimum qualified site depth, up to 4,294,967,295 |
| `--ignore-orphan` | — | Off | Skip paired reads without the SAM proper-pair flag; retain single-end reads |
| `--compression-threads` | Nonnegative integer | `0` | Private BGZF workers per output; with compression, defaults to `1` when multiple caller threads are used |

</div>

Every module accepts one biological sample per BAM. In `joint`, these options
apply to both outputs. See [Calling BAM](../input-output/input-data.md#calling-bam) for input
requirements.

## `bsbit call meth` { #call-meth }

Aggregate strand-specific methylation evidence from a bsbit BAM:

```bash
bsbit call meth \
  -i sample.prep.bam \
  -r reference.fa \
  -o sample.CGmap \
  -f cgmap \
  -t 8
```

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-o`,<br>`--output` | Output path | Required | Output path |
| `-f`,<br>`--format` | `cgmap` or `bed` | Required | Output format: CGmap or extended bedMethyl |
| `--cg-only` | — | Off | Omit CHG and CHH sites |

</div>

See [Call methylation](../guides/methylation-calling.md) for filtering behavior and
[CGmap](file-formats.md#cgmap-methylation-calls) or [extended
bedMethyl](file-formats.md#extended-bedmethyl-calls) for output fields.

## `bsbit call snp` { #call-snp }

Call quality-weighted, bisulfite-aware diploid SNVs:

```bash
bsbit call snp \
  -i sample.prep.bam \
  -r reference.fa \
  -o sample.vcf \
  -t 8
```

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-o`,<br>`--output` | VCF path | Required | VCF output path |
| `--min-alt-count` | Positive integer | `2` | Candidate threshold and selected-ALT `LowAD` threshold, up to 4,294,967,295 |
| `--min-alt-fraction` | Decimal from `0` to `1` | `0.1` | Minimum strongest-ALT candidate fraction |
| `--min-gq` | Integer from `0` to `99` | `0` | `LowGQ` filter threshold; 0 disables it |
| `--min-aq` | Integer from `0` to `99` | `30` | Per-ALT posterior-presence `LowAQ` threshold |
| `--heterozygosity` | Decimal greater than `0` and less than `1` | `0.001` | Reference-divergence prior |
| `--underconversion-rate` | Decimal from `0` to `1` | `0.0025` | Non-conversion probability |
| `--overconversion-rate` | Decimal from `0` to `1` | `0` | Overconversion probability |
| `--sample-name` | Sample name | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |

</div>

See [Call SNVs](../guides/variant-calling.md) for BAM preparation,
the likelihood model, filters, and VCF fields.

## `bsbit call joint` { #call-joint }

Produce methylation and VCF outputs while sharing the first evidence pass.
This command accepts all shared calling options listed above, including base
quality, MAPQ, minimum depth, regions, compression, and worker controls. The
table below lists the additional joint-calling options:

```bash
bsbit call joint \
  -i sample.prep.bam \
  -r reference.fa \
  -p sample \
  -f cgmap \
  -t 8
```

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-p`,<br>`--prefix` | Output prefix | Required | Prefix for the methylation and VCF output paths |
| `-f`,<br>`--meth-format` | `cgmap` or `bed` | Required | Methylation output format: CGmap or extended bedMethyl |
| `--cg-only` | — | Off | Omit CHG and CHH sites from the methylation output |
| `--min-alt-count` | Positive integer | `2` | SNV candidate threshold and selected-ALT `LowAD` threshold, up to 4,294,967,295 |
| `--min-alt-fraction` | Decimal from `0` to `1` | `0.1` | Minimum strongest-ALT candidate fraction |
| `--min-gq` | Integer from `0` to `99` | `0` | `LowGQ` threshold; zero disables the filter |
| `--min-aq` | Integer from `0` to `99` | `30` | Per-ALT posterior-presence `LowAQ` threshold |
| `--heterozygosity` | Decimal greater than `0` and less than `1` | `0.001` | Reference-divergence prior |
| `--underconversion-rate` | Decimal from `0` to `1` | `0.0025` | Non-conversion probability |
| `--overconversion-rate` | Decimal from `0` to `1` | `0` | Overconversion probability |
| `--sample-name` | Sample name | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |

</div>

Base-quality, MAPQ, minimum-depth, region, compression, and orphan settings
apply to both outputs. `--cg-only` affects only methylation; SNV-specific
options affect only variants. With prefix `sample`, the outputs are
`sample.CGmap` (or `sample.bed`) and `sample.vcf`; `-c true` appends `.gz` to
both. `--compression-threads` is per output. See [Call SNVs](../guides/variant-calling.md#run-joint-calling).

## `bsbit combine`

Join named CGmap and/or extended bedMethyl samples into a BED6-plus-matrix
output:

```bash
bsbit combine \
  -i tumor.CGmap,normal.CGmap \
  --sample-name tumor,normal \
  -o cohort.bed \
  --matrix both \
  --min-count 10 \
  --min-prop 0.8 \
  -t 8
```

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | Comma-separated file paths | Required once | One list of sorted 8-column CGmap or 18-column bsbit extended bedMethyl inputs; formats may be mixed |
| `-o`,<br>`--output` | Output path | Required | Exact destination, or filename template when `--matrix both` |
| `-c`,<br>`--compress` | `true` or `false` | `false` | Write deterministic BGZF when `true`; otherwise write plain text |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of hierarchical input-merge workers |
| `--sample-name` | Comma-separated sample names | Input paths | Optional sample labels, one per input |
| `--matrix` | `level`, `count`, or `both` | `level` | Level matrix, count matrix, or both as separate files from one merge |
| `--min-count` | Positive integer | `1` | Per-sample minimum methylated-plus-unmethylated coverage |
| `--min-prop` | Decimal from `0` to `1` | `0` | Minimum proportion of samples passing `--min-count`; at least one is always required |
| `--cg-only` | — | Off | Retain only CpG sites |
| `--compression-threads` | Nonnegative integer | `0` | Private BGZF workers per output; with compression, defaults to `1` when multiple merge threads are used |

</div>

With `--matrix both`, bsbit inserts `.level` and `.count` before the output
format suffix. `.level` identifies the methylation-level matrix, while `.count`
identifies the methylated/total-count matrix. See [Build methylation
matrix](../guides/methylation-matrices.md) for input ordering, filtering, and
matrix schemas.

## Exit codes and output behavior

Successful commands exit with code 0. Invalid command-line usage and
unsupported modes exit with code 2; operational failures exit with code 1.

Help, version information, and textual reports are written to standard output.
Indexing and alignment also report progress there. With `align --metrics`,
standard output contains the machine-readable profiling report instead of
progress messages.
