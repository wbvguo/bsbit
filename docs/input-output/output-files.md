# Output files

`bsbit index`, `bsbit align`, and `bsbit call` write directly to their requested
output paths. An existing regular file will be overwritten. If a command fails,
its output may be empty or partial. See [Filesystem and output](../help/troubleshoot.md#output-cannot-be-opened-for-writing) for permission checks and failure behavior.

## Outputs by stage

| Stage | Example output | What to know |
| --- | --- | --- |
| [Index](../guides/indexing.md) | `reference.bsbit` | Reference index used by `bsbit align` |
| [Alignment](#alignment-bam) | `sample.bam` | Input-order BAM; not yet coordinate-sorted or indexable |
| [BAM preparation](../guides/prepare-bam.md) | `sample.prep.bam` + `.bai` | Coordinate-sorted with a study-appropriate duplicate policy |
| [Methylation calling](../guides/methylation-calling.md) | `methylation.bed` or `methylation.CGmap` | Per-site methylation calls |
| [SNP calling](../guides/variant-calling.md) | `variants.vcf` | Variant calls in VCF format |
| [Joint calling](../guides/variant-calling.md) | `sample.CGmap` + `sample.vcf` | One prefix derives both methylation and variant output names |
| [Matrix aggregation](../guides/methylation-matrices.md) | `cohort.level.bed` and/or `cohort.count.bed` | Level or count matrices from sorted methylation calls |

Name-sorted, fixmate, and position-sorted BAM files are intermediate files.
Retain the final prepared BAM and its index, the authoritative FASTA, the bsbit
reference index, and analysis outputs needed for reproducibility.

## Alignment BAM

`bsbit align` writes a BAM file in FASTQ input order. It preserves the
alignment, reference identity, complete read sequence and qualities, and
bisulfite strand information required by bsbit callers. See
[File formats](../reference/file-formats.md#bam-alignments-and-index) for the
SAM fields, tags, and provenance record.

### Records and ordering

By default, alignment writes one primary record per input read. Accepted
placements are mapped, ambiguous results may retain a deterministic low-MAPQ
representative, and reads without a placement are written as unmapped.
Paired-end input produces one record per mate.

`--mapped-only` removes records without an accepted placement but keeps mapped
MAPQ-0 representatives.

### Output contracts

Use `--output-contract minimal` (the default) to write `NM` and `XG`, or
`--output-contract bismark` to also write `MD`, `XM`, and `XR` following the
Bismark-compatible tag schema. The selected contract does not change the
alignment. Sorting and duplicate handling must preserve the standard `@SQ M5`
reference checksums and mapped-record `XG` tags. The informational `@HD bs`
field and standard `@PG` processing history should also be retained.

See [Prepare BAM file](../guides/prepare-bam.md) for sorting and indexing,
[Validate the BAM](../guides/alignment.md#validate-the-bam) for structural
checks, and [Input data](input-data.md#calling-bam)
for calling requirements.

### Alignment metrics

`bsbit align --metrics` writes an optional profiling TSV to standard output.
It is a diagnostic, not a normal workflow result.

## Compression and output

Alignment BAM is always BGZF-compressed. Calling and matrix outputs are plain
text by default; use `-c true` to write BGZF-compressed output. BGZF-compressed
VCF output can be indexed with
[`bcftools index`](https://samtools.github.io/bcftools/bcftools.html#index);
BED-family outputs can be indexed with
[`tabix`](https://www.htslib.org/doc/tabix.html).
