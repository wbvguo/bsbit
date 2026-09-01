# Call SNVs

`bsbit call snp` writes diploid SNVs. `bsbit call joint` writes SNVs and
methylation from one shared evidence pass. Both are deterministic technical
workflows; neither is a substitute for study-specific scientific validation.

!!! warning "Validate calls for your study"
    Benchmark variants against matched WGS, an independent truth set, or an
    established bisulfite-aware caller before production scientific use.

## Before you run

Prepare `sample.analysis.bam`, then verify the shared [calling input
contract](../reference/input-data.md#calling-bam-and-indexed-reference). Use the
same indexed FASTA assembly as alignment. One BAM represents one biological
sample; `--sample-name` changes only the VCF sample label.

Without a region option, callers process the whole BAM dictionary. Repeatable
`--region` values and `--regions-file` are merged into one interval union, so
overlapping targets cannot duplicate a call. See the [CLI
reference](../reference/cli.md#call-snp) for coordinate syntax and defaults.

## Run

### Run SNV calling

```bash
bsbit call snp \
  --input sample.analysis.bam \
  --reference reference.fa \
  --output sample.vcf.gz \
  --compress true \
  --threads 8
```

[`snp` parameters and defaults](../reference/cli.md#call-snp)

### Run joint calling

```bash
bsbit call joint \
  --input sample.analysis.bam \
  --reference reference.fa \
  --meth sample.cgmap.gz \
  --meth-format cgmap \
  --vcf sample.vcf.gz \
  --compress true \
  --threads 8
```

`joint` applies the same base-quality, MAPQ, region, and compression settings
to both outputs. The methylation and VCF destinations must be different.

[`joint` parameters and defaults](../reference/cli.md#call-joint)

## Interpret the output

The caller follows four operational rules:

1. It excludes top-strand C→T and bottom-strand G→A changes from candidate
   discovery, then applies depth, strongest-ALT count, and ALT-fraction gates.
2. It evaluates all ten unordered diploid genotypes using strand, base quality,
   MAPQ, conversion rates, and adaptive methylation integration for C/G.
3. A reference-divergence prior selects the ALT set and informs `QUAL`/`AQ`;
   `GT` dosage and `PL` remain conditional, prior-free comparisons.
4. Overlapping mates contribute once per fragment/site, and completed regions
   stream in deterministic coordinate order without loading all evidence into
   memory.

The [scientific behavior contract](../scientific-contract.md#downstream-calling-semantics)
defines the exact evidence, overlap, and likelihood semantics.

`AQ` and `GQ` answer different questions. `AQ` measures whether a selected ALT
occurs at least once; `GQ` measures confidence in the complete diploid dosage.
A well-supported ALT can therefore have high AQ while heterozygous versus
homozygous-ALT dosage remains uncertain. `AF` is the conditional expected ALT
fraction after selecting the ALT set.

The configured heterozygosity is a reference-divergence prior, not a population
allele frequency. It affects site and ALT selection but does not bias the
conditional `GT` dosage comparison. Model and integration settings are written
to the VCF header.

The module writes VCF 4.3 SNV records:

| Field | Location | Meaning |
|---|---|---|
| `QUAL` | Fixed column | Confidence that the site is not homozygous reference |
| `GT` | FORMAT | Maximum-likelihood dosage for the selected ALT set |
| `GQ` | FORMAT | Confidence in the complete diploid genotype |
| `AQ` | FORMAT | Confidence that each selected ALT is present |
| `PL` | FORMAT | Prior-free genotype likelihoods |
| `DP` | INFO/FORMAT | Quality-filtered depth |
| `AD` | FORMAT | Raw allele depths |
| `IAD` | FORMAT | Bisulfite-informative allele depths |
| `AF` | INFO | Conditional expected ALT fraction |
| `BSI` | INFO | Strand-discrimination summary |
| `BS8` | INFO | Filtered top-strand A,C,G,T counts, then bottom-strand A,C,G,T counts |

An ALT below `--min-alt-count` receives `LowAD`; a genotype below
`--min-gq` receives `LowGQ`; an ALT below `--min-aq` receives `LowAQ`.
By default, PASS requires `AQ >= 30` and has no GQ floor. Add `--min-gq 20`
when PASS must also require a confident complete genotype.

Calling retains bounded regional state rather than a whole-genome candidate
map. Region and likelihood batches adapt to worker count and candidate density;
output stays coordinate-deterministic across `--threads` values.

## Limits and publication

The caller emits biallelic or multiallelic SNVs only. It does not call indels,
assemble haplotypes, use population or known-sites priors, learn conversion
rates, recalibrate qualities, or claim clinical validation.

`joint` shares overlap-collapsed evidence for methylation aggregation and SNV
candidate discovery, then performs a candidate-only likelihood pass. It is an
efficiency and consistency surface, not a different SNV model.

Destinations are create-only. Joint mode publishes both outputs or neither.
`--compress true` writes deterministic BGZF; index the VCF with:

```bash
tabix -p vcf sample.vcf.gz
```

## Next

- Review [`snp`](../reference/cli.md#call-snp) or
  [`joint`](../reference/cli.md#call-joint) parameters
- [Interpret methylation output](methylation.md)
- [Check supported uses and limitations](../known-differences.md)
