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
  -o sample.cgmap.gz \
  -f cgmap \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | `PATH` | Required | Coordinate-sorted and indexed BAM |
| `-r`,<br>`--reference` | `FASTA` | Required | Reference FASTA used to build the alignment index |
| `-o`,<br>`--output` | `PATH` | Required | Path for methylation calls |
| `-f`,<br>`--format` | `cgmap` or `bed` | Required | Methylation output format |
| `-c`,<br>`--compress` | `BOOL` | `true` | Write BGZF-compressed output |
| `-t`,<br>`--threads` | `N` | `1` | Number of regional calling workers |

</div>

## Advanced parameters

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `--region` | `CONTIG:START-END` | All contigs | Genomic region to call |
| `--regions-file` | `BED` | None | BED file containing regions to call; conflicts with `--region` |
| `--min-bq` | `N` | `20` | Minimum base quality for an observation to be counted |
| `--min-mapq` | `N` | `20` | Minimum alignment MAPQ for an observation to be counted |
| `--min-depth` | `N` | `10` | Minimum qualified depth required at a site |
| `--cg-only` | — | Off | Restrict output to CpG sites |
| `--ignore-orphan` | — | Off | Skip paired reads without the SAM proper-pair flag |

</div>

??? note "Target regions"

    A region uses `CONTIG:START-END` with 1-based inclusive coordinates.
    Separate multiple regions with commas:

    ```bash
    --region chr1:1-100000,chr2:200001-300000
    ```

    For many regions, use `--regions-file`, which accepts plain, gzip, or
    BGZF-compressed BED with 0-based half-open coordinates. It cannot be used
    together with `--region`.

??? note "Quality filters"

    Only bases that meet the base-quality threshold and reads whose alignments
    meet the mapping-quality threshold are included in calling.
    `--ignore-orphan` skips paired reads without the SAM proper-pair flag;
    single-end reads are retained. Sites with fewer than `--min-depth`
    qualified observations are omitted. Add `--cg-only` to omit CHG and CHH
    sites.

See the [CLI reference](../reference/cli.md#call-meth) for accepted ranges and
complete parameter details.

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
