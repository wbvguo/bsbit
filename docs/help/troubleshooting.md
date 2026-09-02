# Troubleshoot

bsbit rejects inputs and stops publication when validation, reference identity,
resource, or output requirements are not satisfied. Start with the exact error
on standard error and the matching section below.

## A source build cannot find HTSlib, htscodecs, or libsais

Initialize both submodule levels and rebuild:

```bash
git submodule update --init --recursive
cargo build --locked --release -p bsbit-cli --bin bsbit
```

The repository owns pinned native dependencies under `external/`; installing a
different system HTSlib does not repair a missing submodule checkout.

## The binary exits with an illegal instruction

The tested production target requires a 64-bit Intel or AMD CPU with
x86-64-v3 support, including AVX2 and POPCNT. Check
the flags visible inside the actual VM, container, or WSL environment:

```bash
lscpu | grep -E 'Architecture|Flags'
```

Move the run to a supported host if the required flags are absent or hidden.

## Publication fails under `/mnt/c` on WSL2

Move the reference build and BAM output to the Linux filesystem, for example:

```bash
mkdir -p ~/work/bsbit
```

Large index and output publication intentionally relies on Linux descriptor and
hard-link semantics. Windows-drive 9p mounts can fail this contract even when
ordinary file creation appears to work.

## An output cannot be replaced

bsbit can replace an existing regular-file result, but it rejects a destination
that is a directory or special file. Check the destination type and parent
directory permissions. If the command fails before publication, the previous
result remains unchanged.

## A `.staging` file remains after interruption

The final BAM is published only after the writer completes. If the process was
killed abruptly, a private staging file may remain while the final BAM is
absent. First confirm that no bsbit process still owns it. Preserve it while
investigating the failure; it is not a committed output and should not be fed
to downstream analysis as the result.

## Paired FASTQ names or counts are inconsistent

R1 and R2 must end together and match by ordinal. Names must be identical or
use corresponding `/1` and `/2` suffixes. Common causes are independently
filtered mates, reordering, a truncated gzip file, or incompatible read-name
rewrites.

Validate compression and record counts before rerunning alignment:

```bash
gzip -t sample_R1.fastq.gz
gzip -t sample_R2.fastq.gz
zcat sample_R1.fastq.gz | awk 'END { print NR / 4 }'
zcat sample_R2.fastq.gz | awk 'END { print NR / 4 }'
```

## An index is rejected as corrupt or stale

Preserve the exact error for diagnosis. Rebuild the index from the trusted
original FASTA with `bsbit index`; only after that single command succeeds,
rerun `bsbit align`. Alignment deliberately never repairs missing, corrupt, or
mismatched index data.

## A reference FASTA is rejected as ordinary gzip

Reference FASTA supports plain or BGZF-compressed input. A `.gz` suffix does
not identify which gzip variant was used; bsbit checks the file content and
rejects ordinary gzip because it cannot provide random access.

Recompress the file as BGZF and create the sidecar indexes required by calling:

```bash
gzip -cd reference.fa.gz | bgzip -c > reference.bgzf.fa.gz
samtools faidx reference.bgzf.fa.gz
```

This creates `reference.bgzf.fa.gz.fai` and
`reference.bgzf.fa.gz.gzi`. Alternatively, use an uncompressed FASTA; its
`.fai` is optional because bsbit can build an in-memory position table.

## Alignment produces fewer mapped records than expected

Ambiguous or unmapped reads do not count as successful mappings. Inspect the
BAM flags and mapping summary rather than comparing only the number of mapped
rows:

```bash
samtools flagstat sample.bam
```

## A downstream tool cannot index the BAM

bsbit output is input-order, not coordinate-sorted. Sort before indexing:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

## A caller rejects the BAM or reference

Every caller run, including one using `-t 1`, requires a coordinate-sorted BAM
with an adjacent BAI or CSI. The BAM must retain the canonical bsbit `@PG`
record. Prepare paired output as shown in the [BAM-preparation
guide](../guides/prepare-bam.md), then confirm both files:

```bash
samtools quickcheck -v sample.prep.bam
samtools index sample.prep.bam
samtools view -H sample.prep.bam | grep '^@PG'
```

### A single-end BAM is rejected

Current `bsbit align` output with read 1 only uses numeric MAPQ and declares
`caller-compatible-directional-single` or
`caller-compatible-nondirectional-single`, according to the selected library
profile. If a single-end BAM is rejected, inspect `@PG`: an older
`standard-directional-single` BAM predates numeric MAPQ and must be realigned
rather than having its header or scores rewritten.

All three call modules also require a reference genome FASTA. When a
plain FASTA has an adjacent FAI, bsbit uses it. Without FAI, bsbit scans the
plain FASTA once and builds an in-memory position table without writing a
sidecar. BGZF FASTA requires both FAI and GZI; ordinary gzip FASTA is rejected.

Every BAM dictionary contig must have one uniquely named FASTA entry with the
same length; FASTA order may differ. The caller normalizes and hashes those
sequences in BAM order and compares the result with structured BAM provenance,
so a same-name, same-length wrong assembly is rejected. It still ignores `MD`.

## A run is slower than the published result

Published timing is scoped to the qualified paired-end aligner (historically measured
under the former `bsbit-align` executable name), host,
dataset, options, storage, and measurement boundary. Check that large files are
on local Linux storage, the same mode, template bounds, and worker counts are
selected, compressed input is read directly, and no competing process is
consuming CPU, memory bandwidth, or storage.

Change one option at a time. See
[alignment settings](../guides/alignment.md#configure-alignment)
and the exact [performance protocol](../performance-evidence.md).
