# Workflow

bsbit separates indexing, alignment, BAM preparation, calling, and cohort
aggregation. Stop after the stage required by the analysis.

## Stages

| Stage | Command | Main input | Output | Guide |
|---|---|---|---|---|
| Index | `bsbit index` | Reference FASTA | Opaque `.bsbit` index | [Build index](../guides/indexing.md) |
| Align | `bsbit align` | Index and FASTQ | Input-order BAM | [Align reads](../guides/alignment.md) |
| Prepare | `samtools` | Alignment BAM | Coordinate-sorted, duplicate-handled as appropriate, indexed BAM | [Prepare a BAM](../guides/prepare-bam.md) |
| Call | `bsbit call` | Prepared BAM and the same indexed FASTA | Methylation output and/or VCF | [Methylation](../outputs/methylation.md) · [SNVs](../outputs/variant-calling.md) |
| Combine | `bsbit combine` | Sorted extended bedMethyl files | Level and/or count matrix | [Build matrix](../outputs/methylation-matrices.md) |

Final destinations must be new local paths. Identical immutable inputs,
options, worker count, and binary produce deterministic output within the
documented contract.

## Supported workflows

| Goal | Stages | Requirements |
|---|---|---|
| Alignment only | Index → Align | Single-end or paired-end FASTQ; coordinate sorting is optional unless another tool requires it |
| Methylation or SNV calling | Index → Align → Prepare → Call | Caller-compatible single- or paired-end alignment and the same reference assembly throughout |
| Cohort methylation matrix | All five stages | One sorted extended bedMethyl input and unique sample name per sample |

Alignment support depends on read layout:

| Input | Status |
|---|---|
| Directional paired-end | Qualified path with caller-compatible provenance and published GRCh38 speed, accuracy, and MAPQ evidence |
| Directional single-end | Caller-compatible numeric MAPQ; published 5M-R1 speed and exact/within-5-bp truth evidence |
| Non-directional paired-end | Four-strand behavior and compatibility are tested; directional benchmark results do not apply |
| Non-directional single-end | Four-strand behavior, global cross-pass classification, and caller compatibility are tested; directional single-end benchmark results do not apply |

Directional paired-end WGBS owns the published qualification. Preprocessed
RRBS and targeted reads are accepted when chemistry and orientation match, but
bsbit does not provide assay-specific trimming or interpretation. PBAT, CRAM,
remote input, and object-store input are unsupported.

Methylation, SNP, and joint calling are deterministic technical workflows that
still require study-specific validation. The SNP caller is diploid, does not
call indels or assemble haplotypes, and makes no clinical-validation claim.

## End-to-end example

Exercise all five stages on the bundled synthetic data from the repository
root:

```bash
bash docs/examples/run-end-to-end.sh
```

The script indexes, aligns, prepares the BAM, calls methylation and one
synthetic SNV, and builds one-sample level and count matrices. It is a smoke
test, not a benchmark. See the [outputs overview](../outputs/index.md) for the
artifact set and the Usage guides above for commands to adapt to real data.
