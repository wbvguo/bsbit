# Call SNVs

`bsbit call snp` calls bisulfite-aware diploid SNVs and writes VCF 4.3. It
reports SNVs only, not indels or haplotypes.

!!! warning "Validate calls for your study"
    Benchmark results against an appropriate truth set or an established
    bisulfite-aware caller before production scientific use.

## Inputs

Calling requires a [prepared BAM](prepare-bam.md) and the same reference FASTA
used for alignment. One BAM represents one biological sample.

## Run SNV calling

```bash
bsbit call snp \
  -i sample.prep.bam \
  -r GRCh38.fa \
  -o sample.vcf.gz \
  -c true \
  -t 8
```

## Common options

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | `PATH` | Required | Prepared coordinate-sorted BAM with an adjacent BAI or CSI |
| `-r`,<br>`--reference` | `FASTA` | Required | Same plain or BGZF-compressed reference FASTA used for alignment |
| `-o`,<br>`--output` | `PATH` | Required | Path for the VCF 4.3 output |
| `--sample-name` | `NAME` | BAM `SM`, then BAM filename stem | Name of the single VCF sample |
| `-c`,<br>`--compress` | `true\|false` | `false` | Whether to write BGZF-compressed VCF |
| `-t`,<br>`--threads` | `N` | `1` | Number of regional calling workers, from 1 to 64 |

For region selection, quality and chemistry parameters, and other advanced
options, see the [`bsbit call snp` CLI reference](../reference/cli.md#call-snp).

## Configure calling

**Target regions.** By default, the caller processes every nonempty BAM contig.
Use repeatable `--region CONTIG:START-END` values, `--regions-file BED`, or both
to restrict calling; overlapping targets are merged.

**Candidate and quality filters.** Use `--min-depth`, `--min-alt-count`, and
`--min-alt-fraction` to control candidate discovery. `--min-aq` filters ALT
presence confidence, while `--min-gq` filters complete-genotype confidence.
Keep the defaults unless alternative thresholds have been validated for the
study.

## Run joint calling { #run-joint-calling }

`bsbit call joint` writes methylation calls and SNVs from one shared evidence
pass:

```bash
bsbit call joint \
  -i sample.prep.bam \
  -r GRCh38.fa \
  -m sample.cgmap.gz \
  -f cgmap \
  -v sample.vcf.gz \
  -c true \
  -t 8
```

Shared region, quality, compression, and thread settings apply to both outputs.
The methylation and VCF paths must be different. See the [`bsbit call joint`
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

## Next

- [Call methylation](methylation.md)
- [`bsbit call snp` CLI reference](../reference/cli.md#call-snp)
- [`bsbit call joint` CLI reference](../reference/cli.md#call-joint)
- [VCF format](../reference/file-formats.md#vcf-variant-calls)
- [Known differences and limitations](../known-differences.md)
