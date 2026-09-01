# Troubleshoot problems

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

## The output already exists

Index, SAM, BAM, CGmap, bedMethyl, VCF, and matrix destinations must be
new. Choose another output path or move the prior artifact aside intentionally.
bsbit never truncates, appends to, or silently replaces a result.

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
original FASTA under a new output path with `bsbit index`; only after that
single command succeeds, rerun `bsbit align`. Alignment deliberately never
creates or overwrites missing, corrupt, or mismatched index data.

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
samtools quickcheck -v sample.analysis.bam
samtools index sample.analysis.bam
samtools view -H sample.analysis.bam | grep '^@PG'
```

### A single-end BAM is rejected

Current `bsbit align` output with read 1 only uses numeric MAPQ and declares
`caller-compatible-directional-single` or
`caller-compatible-nondirectional-single`. If a single-end BAM is rejected,
inspect `@PG`: an older `standard-directional-single` BAM predates numeric MAPQ
and must be realigned rather than having its header or scores rewritten.

All three call modules also require the authoritative reference FASTA. A plain
FASTA needs an adjacent FAI; a BGZF FASTA needs both FAI and GZI. Build the
index with `samtools faidx`:

```bash
samtools faidx reference.fa
```

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
[paired-end alignment settings](../guides/alignment.md#choose-an-alignment-mode)
and the exact [performance protocol](../performance-evidence.md).
