# Call methylation

`bsbit call meth` summarizes methylated and unmethylated observations at CG,
CHG, and CHH sites. It writes site-level calls in CGmap or extended bedMethyl
format.

## Inputs

Calling requires a [prepared BAM](prepare-bam.md) and the same reference FASTA
used for alignment.

## Run methylation calling

```bash
bsbit call meth \
  -i sample.prep.bam \
  -r GRCh38.fa \
  -o sample.CGmap \
  -f cgmap \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | BAM path | Required | Coordinate-sorted and indexed BAM |
| `-r`,<br>`--reference` | FASTA path | Required | Reference FASTA used to build the alignment index |
| `-o`,<br>`--output` | Output path | Required | New path for methylation calls |
| `-f`,<br>`--format` | `cgmap` or `bed` | Required | Methylation output format |
| `-c`,<br>`--compress` | `true` or `false` | `false` | Write BGZF-compressed output when `true` |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of regional calling workers |

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
| `--cg-only` | — | Off | Restrict output to CpG sites |
| `--ignore-orphan` | — | Off | Skip paired reads without the SAM proper-pair flag |

</div>

??? note "Target regions"

    Target regions restrict methylation calling to genomic intervals of
    interest. This is useful for capture panels, selected genes or loci.

    Each `--region` selects one interval using 1-based inclusive coordinates;
    repeat the option to select additional intervals:

    ```bash
    --region chr1:1-100000
    ```

    For many intervals, use `--regions-bed`, which accepts plain, gzip, or
    BGZF-compressed BED3+ with 0-based half-open coordinates. Intervals supplied
    by both options are merged before calling, and overlaps do not duplicate
    observations or output calls.

??? note "Quality filters"

    bsbit filters evidence before counting methylated and unmethylated
    observations:

    - **Alignment:** Unmapped, secondary, supplementary, QC-failed, and
      duplicate-marked records never contribute. Alignments with MAPQ 255 or
      MAPQ lower than `--min-mapq` are also excluded. With `--ignore-orphan`,
      paired records lacking the SAM proper-pair flag are skipped.
    - **Base:** A base contributes only when its Phred quality is at least
      `--min-bq`. Where mates overlap, the fragment contributes at most one
      observation at each genomic position.
    - **Site:** A site is written only when at least `--min-depth` qualified
      observations remain. `--cg-only` additionally restricts output to CpG
      sites, omitting CHG and CHH sites regardless of their depth.

See the [CLI reference](../reference/cli.md#call-meth) for accepted ranges and
complete option details.

## Output

Calls follow the BAM contig order and are sorted by genomic position. The two
strands of a CpG are reported separately. In conventional bisulfite sequencing,
an unconverted cytosine is reported as methylated, without distinguishing 5mC
from 5hmC.

See the [CGmap](../reference/file-formats.md#cgmap-methylation-calls) and
[extended bedMethyl](../reference/file-formats.md#extended-bedmethyl-calls)
format descriptions for schemas and coordinate conventions.

## Next

- [Build methylation matrix](methylation-matrices.md)
- [Call SNVs](variant-calling.md)
- [CLI reference](../reference/cli.md#call-meth)
- [File formats](../reference/file-formats.md#cgmap-methylation-calls)
