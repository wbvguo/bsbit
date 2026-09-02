# Quick start

This page provides quick access to commonly used bsbit commands. For input
requirements, complete options, and output details, please check the Usage
guide linked below each command.

## Build a reference index

Build the bsbit alignment index:

```bash
bsbit index \
  -r GRCh38.fa \
  -o GRCh38.bsbit \
  -t 8
```

??? output-behavior "Output behavior"

    ⚠️ Outputs are published only after successful completion. If a file already
    exists at the specified output path, it will be automatically overwritten.

[Usage: Build index](../guides/indexing.md)

## Align reads

For paired-end reads:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample.bam \
  -t 8
```

For single-end reads, supply only read 1:

```bash
bsbit align \
  -i GRCh38.bsbit \
  -1 sample.fastq.gz \
  -o sample.bam \
  -t 8
```

[Usage: Align reads](../guides/alignment.md)

## Prepare BAM file

Before calling, prepare a coordinate-sorted, duplicate-marked, and indexed BAM:

```bash
samtools sort -n -o sample.qname.bam sample.bam
samtools fixmate -m sample.qname.bam sample.fixmate.bam
samtools sort -o sample.sorted.bam sample.fixmate.bam
samtools markdup sample.sorted.bam sample.prep.bam
samtools index sample.prep.bam
```

??? note "When to skip duplicate marking"

    Coordinate-based duplicate marking may remove valid reads from amplicon or
    other fixed-end libraries, where independent molecules can share the same
    coordinates. UMI libraries should use a UMI-aware method instead. When the
    selected duplicate policy does not use `samtools markdup`, coordinate-sort
    and index the BAM directly:

    ```bash
    samtools sort -o sample.prep.bam sample.bam
    samtools index sample.prep.bam
    ```

[Usage: Prepare BAM file](../guides/prepare-bam.md)

## Call methylation or SNVs

For faster repeated calling, optionally build a FASTA index. Calling a plain
FASTA without this file is also supported:

```bash
samtools faidx GRCh38.fa
```

Call methylation:

```bash
bsbit call meth \
  --input sample.prep.bam \
  -r GRCh38.fa \
  --output sample.cgmap.gz \
  --format cgmap \
  --compress true \
  --threads 8
```

[Usage: Call methylation](../guides/methylation.md)

Call SNVs:

```bash
bsbit call snp \
  --input sample.prep.bam \
  -r GRCh38.fa \
  --output sample.vcf.gz \
  --compress true \
  --threads 8
```

[Usage: Call SNVs](../guides/variant-calling.md)

To produce both result types from one evidence pass, use `bsbit call joint`;
the [SNV usage guide](../guides/variant-calling.md#run-joint-calling) has the
command.

## Build methylation matrix

Combine sorted CGmap or extended bedMethyl outputs. Different samples may use
either format:

```bash
bsbit combine \
  --input tumor.cgmap.gz,normal.bed.gz \
  --sample-name tumor,normal \
  --output cohort.bed.gz \
  --matrix both \
  --min-count 10 \
  --min-prop 0.8 \
  --compress true \
  --threads 8
```

[Usage: Build methylation matrix](../guides/methylation-matrices.md)

## Find help

Use `bsbit --help` to list the available top-level commands. Add `--help` after
a command or command module to see the options available at that level:

```bash
bsbit --help
bsbit align --help
bsbit call --help
bsbit call meth --help
```

For help choosing which stages to run, see the [workflow guide](workflow.md).
For a complete option lookup, use the [CLI reference](../reference/cli.md).
