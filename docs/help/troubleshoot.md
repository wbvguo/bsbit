# Troubleshoot

## Installation and platform

### `bsbit: command not found`

The executable is either not installed or its installation directory is not on
`PATH`. If bsbit was installed in a Conda environment, activate that environment
again:

```bash
conda activate bsenv
bsbit --version
```

For a Cargo installation outside Conda, ensure that Cargo's binary directory is
on `PATH`:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
bsbit --version
```

If bsbit was built but not installed, run it from the source tree as
`./target/release/bsbit`, or follow [Build and install
bsbit](../getting-started/installation.md#build-and-install-bsbit).

### A source build cannot find HTSlib, htscodecs, or libsais

Initialize both submodule levels and rebuild:

```bash
git submodule update --init --recursive
cargo build --locked --release -p bsbit-cli --bin bsbit
```

The repository includes pinned native dependencies under `external/`.
Installing a different system HTSlib does not repair missing submodules.

### The binary exits with an illegal instruction

bsbit detects CPU features at runtime and selects the fastest available SIMD
backend. On x86-64, it tries AVX-512, AVX2+POPCNT, SSE4.2+POPCNT, and SSE2 in
that order; the AVX-512 backend requires AVX-512F, AVX-512BW, AVX2, and POPCNT.
AArch64 uses its NEON baseline. A portable scalar backend is available for
diagnosing SIMD problems.

If an illegal-instruction error occurs, first confirm that the executable and
the current system use the same architecture. Then inspect the CPU features
visible to the process, especially inside a VM, container, WSL environment, or
batch allocation:

```bash
lscpu | grep -E 'Architecture|Flags'
bsbit cpu
```

If `bsbit cpu` also exits with an illegal instruction, first confirm which
executable is running:

```bash
command -v bsbit
```

Rebuild it from the current source without `target-cpu=native`, global
`target-feature`, or `-march` flags. Requesting a backend that the current CPU
does not support should return an error. If it causes an illegal instruction
instead, report it as a bug.

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
[Paired-read synchronization](../input-output/input-data.md#paired-read-synchronization).

### An index is rejected as corrupt or stale

Rebuild the index from the trusted reference FASTA, then rerun alignment.
`bsbit align` does not repair or modify an index:

```bash
bsbit index -r reference.fa -o reference.bsbit
```

### A reference FASTA is rejected as ordinary gzip

Index construction can decode ordinary gzip, but `bsbit call` requires
uncompressed or BGZF FASTA for random access. Convert ordinary gzip to BGZF
and create the sidecars required for calling:

```bash
gzip -cd downloaded-reference.gz | bgzip -c > reference.fa.gz
samtools faidx reference.fa.gz
```

This creates `.fai` and `.gzi` files. Alternatively, use an uncompressed
FASTA. See [FASTA reference](../input-output/input-data.md#fasta-reference).

### An uncompressed FASTA reports that `.fai` is missing

The call continues after scanning the uncompressed FASTA once to build a
temporary in-memory line-layout index. Create the adjacent `.fai` to avoid this
scan on later commands:

```bash
samtools faidx reference.fa
```

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
samtools sort -o sample.prep.bam sample.bam
samtools index sample.prep.bam
```

### A caller rejects the BAM or reference

Confirm that:

- the BAM is coordinate-sorted and has a BAI or CSI;
- the standard `@SQ M5` fields and mapped-record `XG` tags remain; and
- the reference is the same FASTA used to build the alignment index.

Validate the BAM and recreate its index after final processing:

```bash
samtools quickcheck -v sample.prep.bam
samtools index sample.prep.bam
samtools view -H sample.prep.bam | grep '^@SQ'
```

For uncompressed FASTA, create the recommended adjacent `.fai`. For BGZF
FASTA, confirm that both `.fai` and `.gzi` exist. The caller also verifies every
standard `@SQ M5` checksum against the FASTA, so matching contig names and
lengths alone is not sufficient. See [FASTA reference](../input-output/input-data.md#fasta-reference).

## Filesystem and output

### Output cannot be opened for writing

This error usually means that the parent directory does not exist, the current
user cannot write there, or the output path is not a regular file. Inspect the
path and its parent:

```bash
output=/path/to/result
output_dir=$(dirname -- "$output")
ls -ld -- "$output_dir"
test ! -e "$output" || ls -l -- "$output"
```

Use an existing directory that you own and choose a normal file path rather
than a directory, symbolic link, or special file. bsbit does not create missing
parent directories, and the output must not refer to one of the command's input
files. Avoid using `sudo`; it can leave root-owned outputs that later runs
cannot overwrite.

### Output under `/mnt/c` on WSL2 is slow or fails

Move the inputs and output to the WSL2 Linux filesystem, such as
`~/work/bsbit`, and rerun the command. Copy completed results back to `/mnt/c`
only when Windows applications need them. `combine` is not qualified for
Windows-mounted paths.

## Performance

### Alignment is slower than expected

Keep large files on local Linux storage and check for competing CPU, memory, or
storage workloads. `--sensitive` and `--non-directional` perform additional
search work, and the mapping-worker count also affects throughput.

Change one setting at a time. See [alignment
settings](../guides/alignment.md#advanced-options) and enable `--metrics`
only when profiling is needed.
