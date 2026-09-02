# CLI reference

This page lists the supported user-facing parameters, aliases, meanings, and
defaults for the `bsbit` umbrella command. `bsbit align` is the standard entry
point for both single-end and paired-end reads. `bsbit` is the only user-facing
executable; there are no separate `bsbit-align` or `bsbit-call` executables.

In the **Default** column, **Required** means that the command has no default,
**Conditional** means that another option or mode determines whether it is
required, and **None** means the option is optional and disabled when omitted.
Existing regular-file destinations are atomically replaced after the new output
is complete. Directories and special files are rejected as destinations.

## Choose a command

| Command | Purpose |
|---|---|
| [`bsbit index`](#bsbit-index) | Build the complete opaque alignment index |
| [`bsbit align`](#bsbit-align) | Standard single-end or paired-end alignment |
| [`bsbit call meth`](#call-meth) | Methylation calling |
| [`bsbit call snp`](#call-snp) | Bisulfite-aware diploid SNV calling |
| [`bsbit call joint`](#call-joint) | Shared methylation and SNV calling |
| [`bsbit combine`](#bsbit-combine) | Combine CGmap and/or extended bedMethyl samples into matrices |

## Help and version parameters

| Option | Value | Default | Description |
|---|---|---|---|
| `-h`,<br>`--help` | — | None | Print the relevant help for any executable or subcommand and exit |
| `help` | — | None | Positional help alias for top-level `bsbit`, `bsbit call`, and `bsbit combine`; prefer `--help` in scripts |
| `-V`,<br>`--version` | — | None | Print the `bsbit` workspace CLI version and exit |

Parameter values use these conventions:

- `N` is a base-10 nonnegative integer unless the row gives a narrower range.
- `P` is a decimal probability from 0 to 1 with at most nine fractional
  digits; exponent notation is not accepted. `--heterozygosity` excludes both
  endpoints.
- `true|false` requires the literal value; flag-only options such as `--help`
  take no value.
- `A|B` lists the only accepted literal choices.

## `bsbit index`

Create the complete alignment index from local plain or BGZF-compressed FASTA.
Ordinary gzip and other compression formats are rejected based on file content.
An existing index bundle is replaced only after the complete new bundle is ready.

```bash
bsbit index -r PATH -o PATH [-t N] \
  [--index-speed balanced|fast]
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-r`,<br>`--reference` | `PATH` | Required | Plain or BGZF-compressed reference genome FASTA |
| `-o`,<br>`--output` | `PATH` | Required | Generated bsbit alignment index |
| `-t`,<br>`--threads` | `N` | `1` | Indexing workers, 1–64 |
| `--index-speed` | `balanced\|fast` | `balanced` | Sparse suffix-array tradeoff: `fast` reduces locate work while increasing index size and mapping RSS |

This single command builds all alignment data. `OUTPUT` is the opaque index
handle passed to `bsbit align`; physical construction and
layout are not part of the user-facing contract. The default `balanced` index
uses stride 16. `fast` uses stride 8 and remains query-compatible with the same
aligner, but it is a distinct on-disk format revision and consumes more storage
and resident memory.

## `bsbit align` { #bsbit-align }

Align one FASTQ, or synchronized paired FASTQ, with the index created by
`bsbit index` and publish BAM in input order. Read 1 alone selects single-end
layout; supplying both read 1 and read 2 selects paired-end layout:

```bash
bsbit align \
  -i reference.bsbit \
  -1 READS_OR_R1.fastq.gz \
  [-2 R2.fastq.gz] \
  -o OUTPUT.bam \
  [OPTIONS]
```

| Shared option | Value | Default | Description |
|---|---|---|---|
| `-i`,<br>`--index` | `PATH` | Required | Opaque index handle created by `bsbit index` |
| `-1`,<br>`--read1` | `PATH` | Required | The single FASTQ, or R1 in a pair |
| `-2`,<br>`--read2` | `PATH` | None | R2; valid only together with read 1 |
| `-o`,<br>`--output` | `PATH` | Required | BAM destination |
| `--sensitive` | — | Off | Audit the default result against the wider bounded candidate frontier |
| `--non-directional` | — | Off | Make one placement decision across all four bisulfite directions |
| `--output-contract` | `minimal\|bismark` | `minimal` | Emit `NM/XG`, or add Bismark-compatible `MD/XM/XR` tags |
| `--mapped-only` | — | Off | Omit primary records without an accepted placement; retained MAPQ-0 placements remain |
| `-t`,<br>`--threads` | `N` | `1` | Mapping workers, 1–64 |
| `--compression-threads` | `N` | `1` | BGZF output workers; use 0 for synchronous compression |
| `--compression-level` | `default\|0..9` | `1` | HTSlib/BGZF compression setting |
| `--metrics` | — | Off | Write a layout-specific profiling TSV to stdout; normal runs keep stdout clean |

Single-end alignment uses the same persisted combined index and bounded
exact-reference verification core as paired alignment. Directional mode searches
OT/OB; `--non-directional` also searches CTOT/CTOB and merges both passes before
classification. Unique placements carry numeric Q10/Q15/Q20/Q30/Q40 evidence
tiers, while tied origins and unmapped records carry Q0. Output declares the
matching `caller-compatible-directional-single` or
`caller-compatible-nondirectional-single` provenance and is accepted by
`bsbit call` after the documented sort, index, reference, and tag checks.
Pair-specific options fail explicitly when only read 1 is supplied.

| Paired-only option | Value | Default | Description |
|---|---|---|---|
| `--total-threads` | `N` | None | Split one 1–64 physical-core budget between mapping and output according to index speed; conflicts with `--threads` and `--compression-threads` |
| `--batch-pairs` | `N` | `16384` | Input pairs per mapping batch |
| `--alignment-queue-batches` | `N` | `2` | Bounded completed-batch queue depth |
| `--min-template-span` | `N` | `0` | Inclusive minimum template span |
| `--max-template-span` | `N` | `1000` | Inclusive maximum template span |

`--metrics` is optional benchmark diagnostics, not an alignment mode. Paired
TSV rows start with `bsbit-alignment-metrics-v2`. Its `fastq_decode_ns` field
contains FASTQ open/decode/close work but excludes bounded-queue blocking, and
`record_worker_total_ns` includes result observation, record composition, and
composer flush. Optional work columns report directional pair passes,
maximal-suffix lanes and rank operations, sampled-SA locate calls/rows/LF/rank,
candidate starts, verified placements, and compatible pair frontiers. These
are profiling counters, not output-quality scores, and are collected only
when `--metrics` is present. Single-end rows start with
`bsbit-single-alignment-metrics-v2` and additionally bind the output contract,
library profile, and complete-versus-mapped-only policy to the reported read
classes, search work, adapter outcomes, and direct-versus-traceback counts.

Default and `--sensitive` are the only public search modes for either layout.
Single-end sensitive preserves the default result as an incumbent, completes
the six-round, 4,096-hit bounded seed frontier as a confidence audit, and does
not invoke pair-only rescue. Its verifier retains d5 for weak, distance-three,
and unmapped incumbents, while a Q20-or-better result at distance two or less
uses the sufficient d3 audit boundary. A different-origin replacement or new
rescue must be unique at Q20 or above; a lower-confidence conflict retains the
incumbent at Q0. Both modes emit
one primary record for each input read or read pair unless either layout uses
`--mapped-only`. The BAM `@PG` line binds the exact reference semantic digest
and alignment mode. The output still needs coordinate sorting, duplicate
handling, and indexing before calling; “caller-compatible” does not mean
“already coordinate-ready.” See the [alignment
guide](../guides/alignment.md#search-sensitivity).

## Calling options shared by `meth`, `snp`, and `joint`

| Option | Value | Default | Description |
|---|---|---|---|
| `-i`,<br>`--input` | `PATH` | Required | Coordinate-sorted, BAI/CSI-indexed caller-compatible bsbit BAM with structured `@PG` provenance |
| `-r`,<br>`--reference` | `FASTA` | Required | Authoritative reference whose normalized content matches the BAM digest; existing FAI is used, otherwise plain FASTA is scanned in memory; BGZF requires FAI/GZI |
| `--region` | `CONTIG:START-END` | Whole dictionary | Repeatable 1-based inclusive target |
| `--regions-file` | `BED` | None | Plain/gzip/BGZF BED3+ targets, 0-based half-open |
| `-c`,<br>`--compress` | `true\|false` | `false` | Write deterministic BGZF; compressed VCF is tabix-compatible |
| `-t`,<br>`--threads` | `N` | `1` | Regional calling workers, 1–64 |
| `--min-base-quality` | `N` | `15` | Minimum observed-base Phred quality, 0–93 |
| `--min-mapq` | `N` | `20` | Minimum mapping quality, 0–254 |

Every module accepts one biological sample per BAM. Region sources form one
merged union; in `joint`, shared options apply to both outputs.

## `bsbit call meth` { #call-meth }

Aggregate strand-specific methylation evidence from a bsbit BAM:

```bash
bsbit call meth \
  --input sample.prep.bam \
  -r reference.fa \
  --output sample.cgmap.gz \
  --format cgmap \
  --compress true \
  --threads 8
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-o`,<br>`--output` | `PATH` | Required | Output path |
| `-f`,<br>`--format` | `cgmap\|bed` | Required | 8-column CGmap or 18-column extended bedMethyl |

All call modules require a coordinate-sorted BAM with BAI/CSI, the matching
authoritative FASTA, and one biological sample per BAM. See [Call
methylation](../guides/methylation.md) for input checks, filtering, and output
semantics.

## `bsbit call snp` { #call-snp }

Call quality-weighted, bisulfite-aware diploid SNVs:

```bash
bsbit call snp \
  --input sample.prep.bam \
  -r reference.fa \
  --output sample.vcf.gz \
  --compress true \
  --threads 8
```

| Option | Value | Default | Description |
|---|---|---|---|
| `--sample-name` | `NAME` | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |
| `-o`,<br>`--output` | `PATH` | Required | VCF path |
| `--min-depth` | `N` | `4` | Minimum candidate and likelihood depth, 1–4,294,967,295 |
| `--min-alt-count` | `N` | `2` | Candidate threshold and selected-ALT `LowAD` threshold, 1–4,294,967,295 |
| `--min-alt-fraction` | `P` | `0.1` | Minimum strongest-ALT candidate fraction, decimal 0–1 |
| `--min-gq` | `N` | `0` | `LowGQ` filter threshold, 0–99; 0 disables it |
| `--min-aq` | `N` | `30` | Per-ALT posterior-presence `LowAQ` threshold, 0–99 |
| `--heterozygosity` | `P` | `0.001` | Reference-divergence prior, decimal strictly between 0 and 1 |
| `--underconversion-rate` | `P` | `0.0025` | Non-conversion probability, decimal 0–1 |
| `--overconversion-rate` | `P` | `0` | Overconversion probability, decimal 0–1 |

See [Call SNVs](../guides/variant-calling.md) for BAM preparation,
the likelihood model, filters, and VCF fields.

## `bsbit call joint` { #call-joint }

Produce methylation and VCF outputs while sharing the first evidence pass:

```bash
bsbit call joint \
  --input sample.prep.bam \
  -r reference.fa \
  --meth sample.cgmap.gz \
  --meth-format cgmap \
  --vcf sample.vcf.gz \
  --compress true \
  --threads 8
```

| Option | Value | Default | Description |
|---|---|---|---|
| `--sample-name` | `NAME` | Unique BAM `SM`, then BAM filename stem | Rename the one VCF sample column |
| `-m`,<br>`--meth` | `PATH` | Required | Methylation-output path |
| `-f`,<br>`--meth-format` | `cgmap\|bed` | Required | 8-column CGmap or 18-column extended bedMethyl |
| `-v`,<br>`--vcf` | `PATH` | Required | VCF 4.3 output path; must differ from the methylation path |
| `--min-depth` | `N` | `4` | Minimum SNV candidate and likelihood depth, 1–4,294,967,295 |
| `--min-alt-count` | `N` | `2` | SNV candidate threshold and selected-ALT `LowAD` threshold, 1–4,294,967,295 |
| `--min-alt-fraction` | `P` | `0.1` | Minimum strongest-ALT candidate fraction, decimal 0–1 |
| `--min-gq` | `N` | `0` | `LowGQ` threshold, 0–99; zero disables the filter |
| `--min-aq` | `N` | `30` | Per-ALT posterior-presence `LowAQ` threshold, 0–99 |
| `--heterozygosity` | `P` | `0.001` | Reference-divergence prior, decimal strictly between 0 and 1 |
| `--underconversion-rate` | `P` | `0.0025` | Non-conversion probability, decimal 0–1 |
| `--overconversion-rate` | `P` | `0` | Overconversion probability, decimal 0–1 |

Base-quality and MAPQ thresholds apply to both evidence streams; the remaining
quality and chemistry options control SNV calling. Regions and compression are
shared by both outputs. The two paths must differ; existing files are replaced
as one unit. See
[Call SNVs](../guides/variant-calling.md) for preprocessing and
output semantics.

## `bsbit combine`

Join named CGmap and/or extended bedMethyl samples into a BED6-plus-matrix
output:

```bash
bsbit combine \
  --input tumor.cgmap.gz,normal.bed.gz \
  --sample-name tumor,normal \
  --output cohort.bed.gz \
  --matrix both \
  --min-count 10 \
  --min-prop 0.8 \
  --compress true \
  --threads 8
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-i`,<br>`--input` | `PATH[,PATH...]` | Required, repeatable | Comma-separated sorted 8-column CGmap or 18-column bsbit extended bedMethyl inputs; formats may be mixed across samples |
| `--sample-name` | `NAME[,NAME...]` | Input paths | Optional comma-separated sample labels, one per input |
| `-o`,<br>`--output` | `PATH` | Required | Destination; with `both`, a template for `.level` and `.count` destinations |
| `-m`,<br>`--matrix` | `level\|count\|both` | `level` | Level matrix, count matrix, or both as separate files from one merge |
| `--min-count` | `N` | `1` | Per-sample minimum methylated-plus-unmethylated coverage |
| `--min-prop` | `P` | `0` | Minimum proportion of samples passing `--min-count`; at least one is always required |
| `-c`,<br>`--compress` | `true\|false` | `false` | Write deterministic BGZF |
| `-t`,<br>`--threads` | `N` | `1` | Hierarchical input merge workers, 1–64 |

With `--matrix both`, the output path is a template for separate `.level` and
`.count` files; the unsuffixed path is not created. `.level` or `.count` is
inserted before `.bed.gz`, `.bed.bgz`, `.bed`, `.gz`, or `.bgz`; for any other
name it is appended. See [Build methylation matrix](../guides/methylation-matrices.md)
for schemas, missing values, filtering, and resource behavior.

## Exit and publication behavior

The CLI uses exit code 2 for usage or unsupported-mode errors and 1 for
operational failures. A completed result atomically replaces an existing
regular-file destination; a failure preserves the previous result. Paired
outputs are replaced and rolled back together. Standard output is reserved for
help and textual reports; result data is written to the selected file paths.
