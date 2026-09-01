# Outputs overview

bsbit creates a final output only after the command succeeds and never
overwrites an existing destination.

## Outputs by stage

- **Index:** `bsbit index` creates an opaque `.bsbit` index used directly by
  `bsbit align`.
- **Alignment:** `bsbit align` creates an input-order BAM. It is valid BAM but
  is not coordinate-sorted or indexable yet.
- **Metrics:** `bsbit align --metrics` writes a two-row TSV to stdout. Redirect
  it explicitly when needed.
- **BAM preparation:** `samtools` creates a coordinate-sorted,
  duplicate-marked BAM plus BAI or CSI.
- **Methylation calling:** `bsbit call meth` creates CGmap or 18-column
  extended bedMethyl.
- **SNP calling:** `bsbit call snp` creates VCF 4.3.
- **Joint calling:** `bsbit call joint` creates methylation output and VCF from
  one evidence pass; both files publish or neither does.
- **Cohort merge:** `bsbit combine` creates level, count, or paired matrix
  files from sorted extended bedMethyl inputs.

## Files from the end-to-end example

- Reference: `reference.fa`, `reference.fa.fai`, and `reference.bsbit`
- Alignment: `alignment.bam` and `alignment.summary.tsv`
- Prepared BAM: `alignment.analysis.bam` and its `.bai`
- Calls: `methylation.bed` and `variants.vcf`
- Cohort matrices: `cohort.level.bed` and `cohort.count.bed`

Name-sorted, fixmate, and position-sorted BAM files are intermediate files.
Retain the final analysis BAM and index, authoritative FASTA and index, bsbit
reference index, and analysis outputs needed for reproducibility.

## Interpret the outputs

- [Alignment BAM](alignment-bam.md): provenance, record completeness, MAPQ,
  standard fields, and auxiliary tags.
- [Prepare a BAM for calling](../guides/prepare-bam.md): sorting, fixmate,
  duplicate handling, and indexing.
- [Methylation output](methylation.md): CGmap and extended bedMethyl schemas.
- [SNP and joint output](variant-calling.md): VCF fields, filters, and quality.
- [Methylation matrices](methylation-matrices.md): sample columns, missing
  values, filtering, and paired-output naming.

## Compression and publication

Alignment BAM is BGZF-compressed. Calling and matrix commands write plain text
unless `--compress true` is selected; compression is not inferred from the
filename. Compressed VCF and BED-family outputs can be indexed with tabix.

All destinations are create-only. A failed command does not publish a final
target. Joint calling and `combine --matrix both` publish their paired outputs
as one unit.
