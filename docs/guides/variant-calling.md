# Call SNVs

`bsbit call snp` identifies genetic variants from bisulfite sequencing data.
Currently, it reports only diploid single-nucleotide variants in VCF format.

## Inputs

Calling requires a [prepared BAM](prepare-bam.md) and the same reference FASTA
used for alignment. One BAM represents one biological sample.

## Run SNV calling

```bash
bsbit call snp \
  -i sample.prep.bam \
  -r GRCh38.fa \
  -o sample.vcf \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | BAM path | Required | Coordinate-sorted and indexed BAM |
| `-r`,<br>`--reference` | FASTA path | Required | Reference FASTA used to build the alignment index |
| `-o`,<br>`--output` | VCF path | Required | New path for the VCF output |
| `-c`,<br>`--compress` | `true` or `false` | `false` | Write BGZF-compressed output when `true` |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of regional calling workers |
| `--sample-name` | Sample name | BAM `SM` or filename stem | VCF sample name |

</div>

## Advanced options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `--region` | `CONTIG:START-END` | All contigs | One 1-based inclusive genomic region |
| `--regions-bed` | BED path | — | BED3+ file for multiple regions |
| `--min-bq` | Integer from `0` to `93` | `20` | Minimum base quality for an observation to be counted |
| `--min-mapq` | Integer from `0` to `254` | `20` | Minimum alignment MAPQ for an observation to be counted |
| `--min-depth` | Positive integer | `10` | Minimum qualified depth required at a site |
| `--min-alt-count` | Positive integer | `2` | Minimum informative ALT count for candidate selection and `PASS` |
| `--min-alt-fraction` | Decimal from `0` to `1` | `0.1` | Minimum fraction of qualified depth supporting the strongest ALT |
| `--min-aq` | Integer from `0` to `99` | `30` | Minimum ALT-presence quality required for `PASS` |
| `--min-gq` | Integer from `0` to `99` | `0` | Minimum genotype quality required for `PASS`; `0` disables `LowGQ` |
| `--heterozygosity` | Decimal greater than `0` and less than `1` | `0.001` | Prior probability that a site differs from the reference |
| `--underconversion-rate` | Decimal from `0` to `1` | `0.0025` | Probability that an unmethylated cytosine remains unconverted |
| `--overconversion-rate` | Decimal from `0` to `1` | `0` | Probability that a methylated cytosine is converted |
| `--ignore-orphan` | — | Off | Skip paired reads without the SAM proper-pair flag |

</div>

??? note "Target regions"

    Target regions restrict methylation calling to genomic intervals of
    interest. This is useful for capture panels or selected genes.
    Use `--region` with 1-based inclusive coordinates for a single region:

    ```bash
    --region chr1:1-100000
    ```

    For many intervals, use `--regions-bed`, which accepts plain, gzip, or
    BGZF-compressed BED3+ with 0-based half-open coordinates. Intervals supplied
    by both options are merged before calling, and overlaps do not duplicate
    observations or output calls.

??? note "Quality filters"

    Only bases that meet the base-quality threshold and reads whose alignments
    meet the mapping-quality threshold are included in calling.
    `--ignore-orphan` skips paired reads without the SAM proper-pair flag;
    single-end reads are retained.
    `--min-depth` sets the qualified depth required for candidate evaluation,
    while `--min-gq` adds `LowGQ` when genotype confidence is below its
    threshold. Filtered records remain in the VCF.

??? note "Allele thresholds"

    The strongest non-reference allele must meet both `--min-alt-count` and
    `--min-alt-fraction` before a site is evaluated as an SNV candidate.
    After genotyping, each selected ALT must also have at least
    `--min-alt-count` bisulfite-informative observations; otherwise the call is
    marked `LowAD`.

    Allele quality (`AQ`) is the Phred-scaled confidence that a selected ALT is
    present. A selected ALT below `--min-aq` adds `LowAQ`. Increasing these
    thresholds requires stronger ALT evidence but may reduce sensitivity.

??? note "Conversion rates and prior"

    `--underconversion-rate` is the probability that an unmethylated cytosine
    fails to convert and is observed as C instead of T. `--overconversion-rate`
    is the probability that a methylated cytosine converts unexpectedly and is
    observed as T instead of C. The caller uses both rates to distinguish
    conversion errors from SNV evidence.

    `--heterozygosity` sets the prior probability that a site differs from the
    reference. Use validated assay-specific estimates when the default rates
    are not appropriate.

See the [CLI reference](../reference/cli.md#call-snp) for accepted inputs and
complete option details.

## Run joint calling { #run-joint-calling }

`bsbit call joint` writes methylation calls and SNVs from one pass:

```bash
bsbit call joint \
  -i sample.prep.bam \
  -r GRCh38.fa \
  -p sample \
  -f cgmap \
  -t 8
```

Shared region, quality, compression, and thread settings apply to both outputs.
Add `--cg-only` to limit only the methylation output to CpG sites. The
prefix produces `sample.CGmap` (or `sample.bed`) and `sample.vcf`; compression
`-c true` adds `.gz` to both names. See the [`bsbit call joint`
CLI reference](../reference/cli.md#call-joint) for joint-specific options.

## Output

The VCF records `GT`, `GQ`, `AQ`, allele depths, and bisulfite-strand evidence.
`AQ` measures confidence that an ALT is present; `GQ` measures confidence in
the complete diploid genotype. See [VCF variant
calls](../reference/file-formats.md#vcf-variant-calls) for the full fields and
filter definitions.

Compressed VCF is tabix-compatible. Create an index when downstream tools
require random access:

```bash
tabix -p vcf sample.vcf.gz
```

Once indexed, use standard VCF tools such as `bcftools` to filter, query,
normalize, or combine calls. The prepared genotype data can then be integrated
with phenotype or methylation measurements for downstream analyses such as
genome-wide association studies (GWAS), allele-specific methylation analysis,
and methylation quantitative trait locus (mQTL) mapping.

## Next

- [Call methylation](methylation-calling.md)
- [CLI reference: `call snp`](../reference/cli.md#call-snp)
- [CLI reference: `call joint`](../reference/cli.md#call-joint)
- [File formats](../reference/file-formats.md#vcf-variant-calls)
- [Limitations and roadmap](../getting-started/workflow.md#limitations-and-roadmap)
