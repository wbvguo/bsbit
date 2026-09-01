# `bsbit align` performance experiment report

## Executive result

Two source changes are worth retaining: a bounded direct mismatch scan and a
correctness-certified single-end record fast path. On this five-million-read
GRCh38 fixture, their combined single-end improvement is about 2--3% wall time
and 1% user CPU relative to the bounded-only build. The bounded scan is neutral
for single-end and saves roughly 1--3% paired-end user CPU.

The larger opportunities are outside those two code paths:

- A stride-8 suffix-array index reduces 8+2 wall time by 6.5% single-end and
  10.4% paired-end, at a 1.443 GiB on-disk SA cost and about 1.44--1.48 GiB
  additional peak RSS. It needs a real runtime-selectable 8/16 format before it
  is mergeable.
- With the stride-16 index, moving the paired pipeline from 8 mapping + 2 BGZF
  workers to 11 + 3 reduces the stable median from 19.73 s to 15.22 s
  (-22.9%), while spending 6.2% more user CPU.
- Stride 8 and 11+3 do not stack on this machine: stride-8 11+3 has a 15.39 s
  median and much higher writer wait. At 14 physical cores the output side is
  already the limiting stage.
- The attempted packed d=3 batch verifier regressed wall time 3.3% and user CPU
  4--5%; the source change was reverted.
- A specialized paired ungapped evaluator had no repeatable gain and increased
  record-worker CPU in its clean comparison; it was reverted.

All valid variants produced byte-identical BAMs: single-end SHA-256
`622a047c6f1ac4999fcc6e5afe8abc2be0d74d7bcb64acf68d9cdbf57743fc6b`
with 5,000,000 records, and paired-end SHA-256
`2d86a537748ffc7a6836bd988b916f44a1bd79b43b15b6c829195af31da7c331`
with 10,000,000 records. The single BAM contains 4,989,023 mapped, 10,977
unmapped, and 4,883,822 `150M` records; paired has 9,884,356 mapped, 115,644
unmapped, and 9,679,831 `150M` records. `samtools quickcheck` was clean, and
these classifications were identical for every checked variant.

## Method

The baseline is clean commit `183d73e62d8256d27864f1c6c4706859e33c48ab`.
The full fixture is five million 150 bp GRCh38 fragments. Production builds use
`x86-64-v3,+popcnt`, fat LTO, one codegen unit, and no incremental compilation.
The standard 8+2 runs are pinned to physical CPUs `0,2,...,18`. Thread-mix runs
use no more than 14 physical cores. The run harness captures wall/user/system
time, peak RSS, the exact binary and command, built-in paired metrics, BAM
SHA-256, record count, `samtools quickcheck`, and record classification.

Runs were A/B interleaved where possible. The Windows/WSL host occasionally
showed large system-time or scheduler excursions while another build was
active. Those raw runs remain on disk and are explicitly marked invalid rather
than folded into medians. The stable conclusions rely on repeated or sandwich
comparisons and user CPU as a secondary signal.

The originally supplied
`/tmp/bsbit-current-benchmark-20260831/indices/bsbit/current.bsbit` has an older,
incompatible catalog magic. Tests therefore use the compatible stride-16 index
at `/tmp/bsbit-flattened-20260831/indices/bsbit/current.bsbit`.

## Retained source experiments

### Bounded complete-distance scan

`ungapped_distance` formerly constructed a complete `UngappedProfile`,
including prefix arrays, only to ask whether the full-span distance was at most
three. The retained implementation scans the already selected reference span
directly, applies the same strand/orientation conversion semantics, and exits
at the fourth mismatch.

Exhaustive two-base differential tests cover every reference/query base pair,
all four bisulfite strands, thresholds 0--3, empty/too-long input, and reference
boundaries. Full-corpus BAM hashes are identical.

The result is intentionally modest. Early single-end A/B runs were 16.46/16.41
s baseline versus 16.56/16.80 s bounded, so there is no demonstrated
single-end wall benefit. Paired user CPU consistently moved from the mid/high
156--157 s range to roughly 151--155 s (about 1--3%), although wall time is
host-noise limited. The final paired profile attributes 4.30% self time to the
new scan, versus 5.50% for profile construction in the earlier profile. This is
a safe low-risk cleanup, not a major speedup.

### Certified single-end direct record construction

For a selected ungapped placement, the retained fast path builds the all-`M`
BAM record directly only when the existing at-most-two mismatch certificate
proves that canonical traceback cannot choose an equally scoring shifted-gap
path. Distance three, indels, reads above the fixed hot-path length, and
shifted-gap ties fall back to the original canonical traceback. This condition
is important: merely seeing a no-indel candidate is not enough to guarantee an
identical CIGAR.

Tests compare the complete `AlignmentRecord` against the old builder at d=0,
d=1, d=2, and reverse orientation, and prove that an `AC`/`CA` shifted-gap tie
falls back. An exhaustive wrapper test covers all two-base sequences and four
strands. Full-corpus output is byte-identical.

Two clean fast-path runs were 15.91 s / 127.85 user s and 15.90 s / 127.43 user
s. A nearby bounded-only run was 16.37 s / 129.01 user s. This supports a
conservative 2--3% wall and about 1% user-CPU improvement. A second sandwich
pair suffered abnormal WSL system time and is retained as invalid. The earlier
8--10% estimate was too optimistic because the existing code already had a
cheap d<=2 certificate before full traceback.

## Reverted negative experiments

### Specialized paired ungapped evaluator

The experiment replaced paired record inspection's profile/CIGAR evaluator
with a direct conversion-aware classification loop. Exhaustive differential
tests and full BAM hashes passed. The clean direct run was 19.41 s / 152.83
user s against 19.53 s / 151.18 user s for the bounded build, while record
worker CPU increased from 18.27 s to 19.35 s. A second run had abnormal host
system time. There is no repeatable benefit, so the source was reverted.

The final paired profile still spends 4.45% in `evaluate_cigar`, which makes
the site visible but does not invalidate the A/B result: a successful rewrite
must improve data flow or remove work rather than merely spell out the same
per-base loop.

### Packed d=3 verifier loads

The experiment packed the 28 reference codes for each lane into four unaligned
64-bit loads, then used AVX2 shuffle/movemask operations in place of scalar
loads. Differential tests passed and the BAM hash remained identical. Against
the original-kernel sandwich, single-end moved from 16.74 s / 130.67 user s to
17.30 s / 135.28--137.23 user s: +3.3% wall and +4--5% CPU. Paired user CPU
moved from 153.14 s to 158.21 s (+3.3%). The extra packing/shuffle work costs
more than the scalar loads it replaces, so the source was reverted.

The final profile confirms that d=3 verification remains the largest single
single-end leaf (22.61%) and 7.64% paired-end. A future attempt should change
the upstream candidate/reference layout so lanes are already contiguous and
vector-ready; repacking inside the hot kernel cannot pay for itself.

## Stride-8 suffix-array experiment

For experiment isolation, both the builder and reader constants were changed
from 16 to 8, a new full GRCh38 index was built, and a matching binary was used.
Those constants were then restored; this deliberately is not committed as a
format change.

The build took 11:43.17 wall, 2,896.98 user, and 76.24 system seconds with
8,964,080 KiB peak RSS. The stride-8 SA sidecar is 4,068,422,848 bytes versus
2,518,547,488 bytes at stride 16: an increment of 1,549,875,360 bytes, or
**1.443 GiB**. This corrects the earlier estimate: about 2.35 GiB is the entire
existing stride-16 SA file, not the increment. Other index sidecars are
unchanged.

At 8+2, valid single-end medians were 16.03 s / 129.59 user s / about 7.30 GiB
RSS for stride 16 and 14.985 s / 119.83 user s / about 8.74 GiB RSS for stride
8: -6.5% wall, -7.5% CPU, and roughly +1.44 GiB RSS. Paired-end moved from
20.14 s / 158.14 user s / about 7.96 GiB RSS to a 18.045 s / 138.40 user s /
about 9.42 GiB RSS median: -10.4% wall, -12.5% CPU, and roughly +1.46 GiB RSS.
All BAMs are byte-identical.

The paired mapping-worker CPU metric fell from roughly 111.69 s at stride 16 to
90.75--92.07 s at stride 8. Writer-queue wait simultaneously rose from about
0.19 s to 0.80--0.98 s, showing that faster mapping begins exposing output
pressure.

A mergeable design should encode/validate SA stride in a format minor version,
load either stride at runtime, and offer stride 8 as an explicit fast-index
build option. Hard-coding the reader to 8 would silently reject or misread
existing indexes and is not acceptable.

## Paired thread/BGZF scaling

All rows below use the retained source fast paths and the stride-16 index.

| mapping + BGZF | cores | median wall | median user | writer wait | conclusion |
| --- | ---: | ---: | ---: | ---: | --- |
| 8 + 2 | 10 | 19.73 s | 155.70 s | 0.079--0.194 s | reference |
| 10 + 4 | 14 | 18.41 s | 168.13 s | 0.005 s in clean run | unstable 16.85--19.97 s; CPU-expensive |
| 11 + 3 | 14 | 15.22 s | 165.41 s | 0.079--0.081 s | best stable mix, -22.9% wall, +6.2% CPU |
| 12 + 2 | 14 | 18.245 s | 170.02 s | 2.47--2.91 s | two BGZF workers bottleneck |

The stride-8 11+3 combination ran 15.75 and 15.03 s (15.39 s median), versus
15.12 and 15.32 s for stride-16 11+3. Stride 8 still reduced user CPU by about
7%, but writer wait increased to 0.66--0.85 s and wall time did not improve.
Thus 11+3 is the immediate whole-program speed win; stride 8 needs more output
headroom to add to it on this 14-core ceiling.

## Final profiles and next priorities

Representative final-source `perf record -e cpu-clock:u -g` runs captured
129k single-end and 161k paired-end samples with zero lost samples. The major
single-end leaves are d=3 batch verification 22.61%, two-lane SA location
11.61%, FM backward wavefront 8.85%, flexible verification 7.84%, fallback
traceback 5.02%, and the direct certificate 2.62%. The major paired leaves are
two-lane SA location 16.23%, FM backward wavefront 10.96%, d=3 batch
verification 7.64%, LF-row 5.77%, scalar backward extension 5.12%, FASTQ decode
4.78%, CIGAR evaluation 4.45%, and bounded ungapped distance 4.30%.

Recommended order from this evidence:

1. Make 11 mapping + 3 BGZF workers the tuned 14-core configuration (or add an
   automatic physical-core-aware split), with writer-wait telemetry retained.
2. Formalize runtime stride-8/stride-16 compatibility as an optional index
   tradeoff; do not change the hard-coded constant alone.
3. Investigate BAM record composition/BGZF batching, because it caps both 12+2
   and stride-8 11+3.
4. Redesign d=3 candidate/reference data upstream so the hot kernel consumes
   prepacked contiguous lanes. Do not retain the measured in-kernel packer.
5. Keep the bounded scan and certified single fast path, but expect only the
   measured small gains from them.

## Evidence locations

- Selected numbers: `results.tsv`, `paired-metrics.tsv`, `index-results.tsv`
- Reproduction commands: `commands.md`
- All timing, metrics, BAM summaries, and checksums:
  `/tmp/bsbit-align-perf-experiments-20260901/runs`
- Experiment binaries: `/tmp/bsbit-align-perf-experiments-20260901/binaries`
- Stride-8 index and build logs:
  `/tmp/bsbit-align-perf-experiments-20260901/sa8-index`
- Final `perf.data` and text reports:
  `/tmp/bsbit-align-perf-experiments-20260901/perf-final`
