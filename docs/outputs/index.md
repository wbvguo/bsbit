# Output files

bsbit writes each result to a temporary file and moves it to the requested
output path only after the command succeeds. If the output file already exists,
bsbit replaces it atomically.

## Outputs by stage

| Stage | Example output | What to know |
|---|---|---|
| [Index](../guides/indexing.md) | `reference.bsbit` | Opaque index used by `bsbit align` |
| [Alignment](#alignment-bam) | `alignment.bam` | Input-order BAM; not yet coordinate-sorted or indexable |
| [BAM preparation](../guides/prepare-bam.md) | `alignment.analysis.bam` + `.bai` | Coordinate-sorted with a study-appropriate duplicate policy |
| [Methylation calling](../guides/methylation.md) | `methylation.bed` or CGmap | Per-site methylation calls |
| [SNP calling](../guides/variant-calling.md) | `variants.vcf` | Variant calls in VCF format |
| [Joint calling](../guides/variant-calling.md) | Methylation output + VCF | Produces both methylation and variant calls |
| [Matrix aggregation](../guides/methylation-matrices.md) | `cohort.level.bed` and/or `cohort.count.bed` | Level or count matrices from sorted methylation calls |

Name-sorted, fixmate, and position-sorted BAM files are intermediate files.
Retain the final prepared BAM and its index, the authoritative FASTA, the bsbit
reference index, and analysis outputs needed for reproducibility.

## Alignment BAM

`bsbit align` writes a SAM/BAM 1.6 file in FASTQ input order. It preserves the
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
`--output-contract bismark` to also write `MD`, `XM`, and `XR`. The selected
contract does not change the alignment. Sorting and duplicate handling must
preserve the bsbit `@PG` header and mapped-record `XG` tags.

See [Prepare BAM file](../guides/prepare-bam.md) for sorting and indexing,
[Validate the BAM](../guides/alignment.md#validate-the-bam) for structural
checks, and [Input data](../reference/input-data.md#calling-bam-and-reference)
for calling requirements.

??? note "Alignment metrics"

    `bsbit align --metrics` writes an optional profiling TSV to standard
    output. It is a diagnostic, not a normal workflow result.

## Compression and output

Alignment BAM and calling and matrix outputs are BGZF-compressed by default.
Use `-c false` to write plain text. BGZF-compressed VCF output can be indexed
with [`bcftools index`](https://samtools.github.io/bcftools/bcftools.html#index);
BED-family outputs can be indexed with
[`tabix`](https://www.htslib.org/doc/tabix.html).
