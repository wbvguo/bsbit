# Troubleshoot

Start with the exact error written to standard error, then use the matching
section below.

## Installation and platform

### A source build cannot find HTSlib, htscodecs, or libsais

Initialize both submodule levels and rebuild:

```bash
git submodule update --init --recursive
cargo build --locked --release -p bsbit-cli --bin bsbit
```

The repository includes pinned native dependencies under `external/`.
Installing a different system HTSlib does not repair missing submodules.

### The binary exits with an illegal instruction

The production build requires a 64-bit Intel or AMD CPU with x86-64-v3 support,
including AVX2 and POPCNT. Check the features visible inside the active VM,
container, or WSL environment:

```bash
lscpu | grep -E 'Architecture|Flags'
```

Use a supported host if these features are absent or hidden.

## Input and reference

### Paired FASTQ names or counts are inconsistent

R1 and R2 must contain the same number of records in the same order. Names must
be identical or use matching `/1` and `/2` suffixes. Independent filtering,
reordering, truncation, and incompatible name rewriting are common causes.

For gzip or BGZF input, check compression and record counts:

```bash
gzip -t sample_R1.fastq.gz
gzip -t sample_R2.fastq.gz
zcat sample_R1.fastq.gz | awk 'END { print NR / 4 }'
zcat sample_R2.fastq.gz | awk 'END { print NR / 4 }'
```

Matching counts do not prove matching names; bsbit validates every pair. See
[Paired-read synchronization](../reference/input-data.md#paired-read-synchronization).

### An index is rejected as corrupt or stale

Rebuild the index from the trusted reference FASTA, then rerun alignment.
`bsbit align` does not repair or modify an index:

```bash
bsbit index -r reference.fa -o reference.bsbit
```

### A reference FASTA is rejected as ordinary gzip

Reference FASTA must be plain or BGZF-compressed. Convert ordinary gzip to
BGZF and create the sidecars required for calling:

```bash
gzip -cd reference.fa.gz | bgzip -c > reference.bgzf.fa.gz
samtools faidx reference.bgzf.fa.gz
```

This creates `.fai` and `.gzi` files. Alternatively, use a plain FASTA. See
[FASTA reference](../reference/input-data.md#fasta-reference).

## BAM and calling

### Alignment produces fewer mapped records than expected

Ambiguous and unmapped reads do not count as successful mappings. Inspect the
BAM summary rather than comparing only mapped rows:

```bash
samtools flagstat sample.bam
```

### A downstream tool cannot index the BAM

The BAM produced by `bsbit align` follows FASTQ input order. Coordinate-sort it
before creating its index:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

### A caller rejects the BAM or reference

Confirm that:

- the BAM is coordinate-sorted and has an adjacent `.bai`;
- the bsbit `@PG` header and mapped-record `XG` tags remain; and
- the reference is the same FASTA used to build the alignment index.

Validate the BAM and recreate its index after final processing:

```bash
samtools quickcheck -v sample.prep.bam
samtools index sample.prep.bam
samtools view -H sample.prep.bam | grep '^@PG'
```

For BGZF FASTA, confirm that both `.fai` and `.gzi` exist. The caller also
verifies the normalized reference digest, so matching contig names and lengths
alone is not sufficient. See [Calling BAM and
reference](../reference/input-data.md#calling-bam-and-reference).

## Filesystem and output

### Writing output fails under `/mnt/c` on WSL2

Build the reference index and write large outputs on the Linux filesystem, for
example under `~/work/bsbit`. Windows-drive 9p mounts may not provide the file
semantics required to finalize these outputs.

### An output cannot be replaced

The destination must be absent or an existing regular file. Directories and
special files are rejected. Check the destination type and parent-directory
permissions.

### A `.staging` file remains after interruption

A staging file may remain if a process is killed while writing output. Confirm
that no bsbit process is still using it. The staging file is not a completed
result and must not be used for downstream analysis.

## Performance

### Alignment is slower than expected

Keep large files on local Linux storage and check for competing CPU, memory, or
storage workloads. `--sensitive` and `--non-directional` perform additional
search work, and the mapping-worker count also affects throughput.

Change one setting at a time. See [alignment
settings](../guides/alignment.md#advanced-parameters) and enable `--metrics`
only when profiling is needed.
