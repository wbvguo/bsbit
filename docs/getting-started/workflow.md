# Workflow

Build the reference index, process each sample through alignment and calling,
then combine the methylation results into a matrix.

<div class="workflow-flow" role="img" aria-label="From left to right: build the index once, repeat alignment, BAM preparation, and calling for every sample, then combine all sample calls">
  <div class="workflow-node">
    <strong>Index</strong>
    <code>bsbit index</code>
  </div>
  <span class="workflow-arrow" aria-hidden="true">→</span>
  <div class="workflow-repeat">
    <span class="workflow-repeat-label">Each sample</span>
    <div class="workflow-row">
      <div class="workflow-node">
        <strong>Align</strong>
        <code>bsbit align</code>
      </div>
      <span class="workflow-arrow" aria-hidden="true">→</span>
      <div class="workflow-node">
        <strong>Prepare BAM</strong>
        <code>samtools</code>
      </div>
      <span class="workflow-arrow" aria-hidden="true">→</span>
      <div class="workflow-node">
        <strong>Call</strong>
        <code>bsbit call</code>
      </div>
    </div>
  </div>
  <span class="workflow-arrow" aria-hidden="true">→</span>
  <div class="workflow-node">
    <strong>Combine</strong>
    <code>bsbit combine</code>
  </div>
</div>

## Stages

| Stage | Command | Main input | Output |
|---|---|---|---|
| [Index](../guides/indexing.md) | `bsbit index` | Reference FASTA | Reusable `.bsbit` alignment index |
| [Align](../guides/alignment.md) | `bsbit align` | Index and FASTQ | Input-order BAM |
| [Prepare BAM](../guides/prepare-bam.md) | `samtools` | Alignment BAM | Coordinate-sorted, duplicate-handled, indexed BAM |
| Call [methylation](../guides/methylation.md) or [SNVs](../guides/variant-calling.md) | `bsbit call` | Prepared BAM and matching reference | Methylation output and/or VCF |
| [Combine](../guides/methylation-matrices.md) | `bsbit combine` | Sorted per-sample methylation call files | Methylation level and/or count matrix |

## Sequencing data support

bsbit supports both directional and non-directional libraries with either
single-end or paired-end data. Directional mode is the default; use
`--non-directional` for non-directional libraries.

Preprocessed RRBS and targeted reads are accepted when chemistry and
orientation match. See the [performance evidence](../performance-evidence.md)
for the validation scope of each alignment mode.

## Limitations and roadmap

The following capabilities are not available in the current release and may be
added in future:

- PBAT and other library protocols, including assay-specific preprocessing and
  interpretation
- Broader variant calling, including indels and haplotype-aware analysis
- Standard-stream (`-`) input and output
- CRAM and additional output formats
- ARM architectures, including Apple Silicon
