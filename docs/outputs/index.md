# Outputs overview

bsbit stages each result privately and publishes it only after the command
succeeds. An existing result file is then replaced atomically.

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
- **Matrix aggregation:** `bsbit combine` creates level, count, or paired
  matrix files from sorted CGmap and/or extended bedMethyl inputs.

## Files from the end-to-end example

- Reference: `reference.fa` and `reference.bsbit`
- Alignment: `alignment.bam` and `alignment.summary.tsv`
- Prepared BAM: `alignment.analysis.bam` and its `.bai`
- Calls: `methylation.bed` and `variants.vcf`
- Cohort matrices: `cohort.level.bed` and `cohort.count.bed`

Name-sorted, fixmate, and position-sorted BAM files are intermediate files.
Retain the final prepared BAM and its index, the authoritative FASTA, the bsbit
reference index, and analysis outputs needed for reproducibility.

## Interpret the outputs

- [Alignment BAM](alignment-bam.md): provenance, record completeness, MAPQ,
  standard fields, and auxiliary tags.
- [Prepare BAM file](../guides/prepare-bam.md): sorting, fixmate,
  duplicate handling, and indexing.
- [Methylation output](../guides/methylation.md): CGmap and extended bedMethyl schemas.
- [SNP and joint output](../guides/variant-calling.md): VCF fields, filters, and quality.
- [Methylation matrices](../guides/methylation-matrices.md): sample columns, missing
  values, filtering, and paired-output naming.

## Compression and publication

Alignment BAM is BGZF-compressed. Calling and matrix commands write plain text
unless `--compress true` is selected; compression is not inferred from the
filename. Compressed VCF and BED-family outputs can be indexed with tabix.

An existing regular-file destination is replaced only after the new output is
complete. A failed command preserves the previous result. Joint calling and
`combine --matrix both` replace their paired outputs as one unit and roll both
back if publication fails.
