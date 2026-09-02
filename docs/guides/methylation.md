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
  -c true \
  -t 8
```

## Common options

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`,<br>`--input` | `PATH` | Required | Prepared coordinate-sorted BAM with an adjacent BAI or CSI |
| `-r`,<br>`--reference` | `FASTA` | Required | Same plain or BGZF-compressed reference FASTA used for alignment |
| `-o`,<br>`--output` | `PATH` | Required | Path for the CGmap or extended bedMethyl output |
| `-f`,<br>`--format` | `cgmap\|bed` | Required | Output format: 8-column CGmap or 18-column extended bedMethyl |
| `-c`,<br>`--compress` | `true\|false` | `false` | Whether to write BGZF-compressed output |
| `-t`,<br>`--threads` | `N` | `1` | Number of regional calling workers, from 1 to 64 |

For region selection, quality thresholds, and other advanced options, see the
[`bsbit call meth` CLI reference](../reference/cli.md#call-meth).

## Configure calling

**Output format.** Use `-f cgmap` for compact 8-column calls or `-f bed` for
18-column extended bedMethyl with explicit intervals, strand, and evidence
categories. Both formats can be used as input to `bsbit combine`.

**Target regions.** By default, the caller processes every nonempty BAM contig.
Use repeatable `--region CONTIG:START-END` values, `--regions-file BED`, or both
to restrict calling; overlapping targets are merged.

**Quality filters.** Observations require base quality 15 and mapping quality 20
by default. Use `--min-base-quality` and `--min-mapq` to change these thresholds.

## Output

Calls are sorted by BAM dictionary, coordinate, and cytosine strand. Opposite
strands of a CpG dyad remain separate records. In conventional bisulfite data,
“methylated” means protected or unconverted cytosine evidence and does not
distinguish 5mC from 5hmC.

See [CGmap](../reference/file-formats.md#cgmap-methylation-calls) and [extended
bedMethyl](../reference/file-formats.md#extended-bedmethyl-calls) for the full
schemas and coordinate conventions.

## Next

- [Build methylation matrix](methylation-matrices.md)
- [Call SNVs](variant-calling.md)
- [`bsbit call meth` CLI reference](../reference/cli.md#call-meth)
- [Methylation output formats](../reference/file-formats.md#cgmap-methylation-calls)
