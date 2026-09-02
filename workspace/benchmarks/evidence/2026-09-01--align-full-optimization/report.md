# `bsbit align` full optimization report

## Executive result

The branch turns every promising direction from the exploratory profile into a
runtime-safe implementation or a measured rejection. It also finds and retains
one additional FASTQ decode optimization after profiling the combined result.

On the five-million-read single-end GRCh38 fixture, the final fast-index source
has a 13.59 s A/B median versus 16.03 s for the earlier retained stride-16
source, about **15.2% faster wall time**. That comparison includes the optional
stride-8 index. The final ASCII-classifier change alone moves its matched
stride-8 median from 14.315 to 13.590 s (-5.1%) and user CPU from 108.76 to
105.44 s (-3.0%).

For paired-end on the same 14 physical cores, the earlier best stride-16 11+3
median was 15.22 s. The clean combined stride-8 10+4 run is 13.31 s (-12.5%);
the later lookup-table A/B medians are 14.13 s control and 13.98 s final while
the host incurred an abnormal 6.5--9.6 s of system CPU. A conservative reading
is **8--12.5% same-core wall improvement** and about 15% lower user CPU. Against
the old 10-core 8+2 median of 19.73 s, the clean 13.31 s run is 32.5% faster,
but that comparison also spends four more physical cores and is not an
algorithm-only claim.

All checked variants retain byte-identical BAM output:

- single-end SHA-256:
  `622a047c6f1ac4999fcc6e5afe8abc2be0d74d7bcb64acf68d9cdbf57743fc6b`
- paired-end SHA-256:
  `2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331`

Every harness run also passed `samtools quickcheck`, record-count, and mapping /
CIGAR classification checks.

## Retained changes

### Versioned optional stride-8 sparse suffix array

`bsbit index --index-speed fast` now emits a format-minor-1 stride-8 sparse SA;
the default `balanced` mode keeps the existing minor-0 stride-16 layout. The
reader validates the minor/stride pair and the exact sample count, then
dispatches once to stride-specialized locate code. Existing stride-16 indexes
remain readable and no per-LF-step runtime stride branch was added.

The isolated stride experiment measured 6.5% single-end and 10.4% paired-end
wall improvement at 8+2. The SA sidecar grows from 2,518,547,488 to
4,068,422,848 bytes: +1,549,875,360 bytes (**1.443 GiB**). Peak mapping RSS
increases by roughly 1.44--1.48 GiB. This is therefore an explicit throughput /
memory tradeoff, not a silent default change.

### Native BAM call buffer

Successful native HTS calls formerly zeroed a 512-byte error buffer on every
record operation even though the text was only read on failure. The retained
wrapper uses uninitialized storage, while the C boundary always initializes at
least a terminating NUL and Rust reads the buffer only on error.

The paired BAM-write stage median moves from 11.754 to 11.181 s, a **4.9%
stage improvement**. Whole-program wall changed by about 3.5%, but host and
mapping variation make the writer-stage metric the defensible attribution.

### Upstream position-major d=3 candidate slab

The successful d=3 redesign interleaves up to four candidate windows while
copying reference bases upstream. Each DP diagonal then consumes one contiguous
32-bit code group before shuffle/movemask classification. This avoids the
failed exploratory design that repacked row-major inputs inside the hot kernel.

Matched stride-16 results:

- single-end stable wall: 15.85 to a 14.705 s median (-7.2%); user-CPU median
  including the wall-invalid control is 128.56 to 117.91 s (-8.3%);
- paired-end wall median: 15.205 to 14.625 s (-3.8%); user CPU: 163.89 to
  157.12 s (-4.1%); mapping-worker CPU: 114.35 to 107.27 s (-6.2%).

A differential test compares interleaved and row-major kernels across read
lengths and candidate counts. Full-corpus hashes remain identical.

### Index-aware paired thread budget

Paired `bsbit align --total-threads 14` now treats 14 as one physical-core
budget and chooses the qualified split from metadata:

| index | mapping + BGZF | reason |
| --- | ---: | --- |
| balanced / stride 16 | 11 + 3 | mapping remains the limiting stage |
| fast / stride 8 | 10 + 4 | faster locate exposes BAM compression |

With stride 8, 11+3 takes 13.91 s and waits 0.704 s on the writer queue. The
10+4 run takes 13.31 s with 0.004 s writer wait; 9+5 regresses to 14.16 s. With
stride 16, 10+4 is 15.07 s and the earlier 11+3 split remains preferable. The
automatic real runs select 11+3 and 10+4 respectively, and both produce the
qualified BAM digest.

### FASTQ ASCII classification table

The final profile showed sequence normalization consuming 3.40% single-end
self time and paired decode running close to the end-to-end critical path. The
old per-base match compiled to an indirect jump table. A 256-byte compile-time
classifier now makes the valid A/C/G/T/N path one lookup, one range check, and
one byte store while preserving the distinction and exact offset for
unsupported IUPAC and invalid bytes.

The exhaustive HTS parser oracle still tests every embeddable byte. Matched
single-end A/B improves wall 5.1% and user CPU 3.0%. Paired decoder medians move
from 12.477 to 12.323 s (-1.2%); whole-program paired wall changes about 1.1%
with neutral user CPU under unusually high host system time, so no larger
paired claim is made.

## Earlier retained fast paths and revised estimates

This branch starts from the prior measured evidence. The bounded ungapped scan
is safe but neutral for single-end and saves roughly 1--3% paired user CPU. The
certified single-end direct BAM materialization contributes about 2--3% wall
time, not the initial 8--10% estimate. A paired direct evaluator did not improve
CPU or wall time and remains reverted. Those measured revisions are preserved
rather than extrapolating from hotspot percentages.

## Rejected experiments

- Repacking d=3 candidate bytes inside the kernel regressed single wall 3.3%
  and CPU 4--5%; it remains reverted.
- A specialized paired ungapped/CIGAR evaluator had no repeatable gain and
  increased record-worker CPU; it remains reverted.
- Moving the d=3 substitution-mask table from per-call stack construction to a
  static table gives 13.815 s control versus 13.830 s experimental wall
  medians (+0.1%), with only -0.4% user CPU. It was reverted as noise-level.
- Increasing BGZF workers without matching index/mapping pressure either wastes
  mapping cores or exposes a different bottleneck. The retained automatic split
  is metadata-aware rather than globally hard-coded.

## Final profile and practical ceiling

The final frame-pointer profiles use the fast index and capture 104K
single-end plus 145K paired-end `cpu-clock:u` samples, both with zero lost
samples. The single-end sequence normalizer has disappeared from the top 0.25%
leaves. Remaining major self-time is:

| hotspot | single | paired |
| --- | ---: | ---: |
| flexible candidate verification wrapper | 14.56% | 2.99% |
| interleaved d=3 kernel | 11.31% | 3.66% |
| FM backward wavefront | 11.31% | 12.96% |
| two-lane SA locate | 8.89% | 11.30% |
| scalar backward extension | below top group | 6.41% |
| paired FASTQ decode wrapper | n/a | 5.62% |
| paired CIGAR evaluation | n/a | 4.95% |

This is not an absolute performance endpoint, but it is the end of the current
low-risk/local pass. The next material gains require one of three larger
projects: a different FM rank/lookup representation, a deeper SIMD/algorithmic
d=3 verifier redesign, or a columnar/batched FASTQ and BAM record pipeline.
Each changes a core data contract and needs its own profiling branch and
correctness campaign. The visible paired CIGAR leaf alone is not sufficient:
the already measured direct rewrite did not pay.

## Evidence locations

- Selected timing rows: `results.tsv`
- Paired internal metrics: `paired-metrics.tsv`
- Reproduction commands: `commands.md`
- All raw runs and BAMs:
  `/tmp/bsbit-align-full-optimization-20260901/runs`
- Final optimized profiles:
  `/tmp/bsbit-align-full-optimization-20260901/perf-final-optimized`
- Release/profile binaries and formal SA8 fixture:
  `/tmp/bsbit-align-full-optimization-20260901`
