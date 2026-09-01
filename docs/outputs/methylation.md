# Call methylation

`bsbit call meth` counts strand-specific methylated and unmethylated
bisulfite observations at CG, CHG, and CHH sites. It writes CGmap or extended
bedMethyl; it does not infer a corrected biological methylation state.

## Before you run

Create `sample.analysis.bam` with the [BAM preparation
recipe](../guides/prepare-bam.md), then verify the shared [calling input
contract](../reference/input-data.md#calling-bam-and-indexed-reference). Use the
same indexed FASTA assembly as alignment.

## Run

```bash
bsbit call meth \
  --input sample.analysis.bam \
  --reference reference.fa \
  --output sample.cgmap.gz \
  --format cgmap \
  --compress true \
  --threads 8
```

Use `--format bed` for extended bedMethyl. Compression is controlled only by
`--compress`; it is not inferred from the filename. Thread counts from 1 to 64
change regional parallelism but not output order.

[`bsbit call meth` parameters and defaults](../reference/cli.md#call-meth)

Without a region option, every nonempty BAM contig is called. Repeatable
`--region` targets and a `--regions-file` form one merged union, so overlapping
targets cannot be counted twice. Exact coordinate conventions are in the [CLI
reference](../reference/cli.md#call-meth).

## Interpret the output

### CGmap

CGmap is tab-delimited and sorted by BAM dictionary, coordinate, and strand:

| Column | Meaning |
|---:|---|
| 1 | Contig |
| 2 | Forward-reference base, `C` or `G` |
| 3 | 1-based position |
| 4 | `CG`, `CHG`, or `CHH` context |
| 5 | Cytosine-strand dinucleotide |
| 6 | Methylated fraction, rendered to six decimals |
| 7 | Methylated observations |
| 8 | Valid methylated + unmethylated coverage |

At a forward-reference C, C is methylated evidence and T is unmethylated
evidence. At a forward-reference G, G is methylated evidence and A is
unmethylated evidence. Other bases do not enter valid coverage.

### Extended bedMethyl

BED output uses 0-based half-open one-base intervals and 18 columns:

| Columns | Meaning |
|---:|---|
| 1–3 | Contig, start, end |
| 4–6 | Modification/context, valid coverage, strand |
| 7–9 | Display interval and color |
| 10–13 | Coverage, percent, methylated, unmethylated |
| 14–18 | Other modification, deletion, failed, different-base, no-call counts |

The caller reports zero for other modification, failed, and no-call counts.
Deletions and non-C/T or non-G/A observations remain in dedicated columns but
do not enter valid coverage.

Opposite strands of one CpG dyad remain separate coordinates. Aggregate them
downstream only when the study design requires a dyad-level estimate.

## Limits and publication

“Methylated” means protected or unconverted cytosine evidence, not direct proof
of 5mC. Values are not corrected for incomplete conversion, and conventional
bisulfite sequencing cannot distinguish 5mC from 5hmC. Context at contig ends
or near noncanonical bases may remain unresolved and is reported in warnings.

Unmapped, secondary, supplementary, QC-failed, and duplicate-marked records
do not contribute. Overlapping mates contribute at most once per fragment and
site under the [scientific evidence rules](../scientific-contract.md#downstream-calling-semantics).

Output publication is atomic and create-only. Plain and BGZF output is
reproducible for the same inputs, options, and binary.

## Next

- [Build matrix](methylation-matrices.md)
- [Review all methylation parameters](../reference/cli.md#call-meth)
- [Check supported uses and limitations](../known-differences.md)
