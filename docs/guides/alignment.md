# Align reads

`bsbit align` is the entry point for both read layouts. Supplying read 1 alone
selects single-end alignment; adding synchronized read 2 selects paired-end
alignment without changing the command. Directional libraries are the default
for either layout, and `--non-directional` enables four-strand placement.

Both layouts open the opaque index created by `bsbit index`; alignment never
builds or modifies index data.

## Before you run

- Complete the [installation](../getting-started/installation.md) and make the
  built binary available on `PATH`.
- [Build one complete index](indexing.md) from the reference FASTA.
- For paired input, confirm that R1 and R2 end together and use synchronized
  read names.
- Choose template-span bounds appropriate to a paired library rather than
  relying on a learned insert-size distribution.

The audited fat-LTO build used for frozen benchmark reproduction is documented
with the [performance evidence](../performance-evidence.md#reproduction-protocol).
It is not required for an ordinary source build or analysis run.

## Align single-end reads

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample.fastq.gz \
  --output-bam sample.bam \
  --sensitive \
  --threads 8
```

Omit `--sensitive` for the low-latency default. With it, single-end alignment
first retains the default result as an incumbent, replays seed intervals
admitted only by the wider 4,096-hit bound, completes the six-round bounded
frontier, and performs distance-5 verification. A different-origin replacement
or new rescue must be unique at MAPQ 20 or above. A lower-confidence conflict
keeps the incumbent as ambiguous at MAPQ 0, and an uncertified rescue remains
unmapped. Pair geometry, mate rescue, and paired adapter recovery do not apply
to a single read.

Directional default and sensitive output also inspect the final 30-base domain
for an exact supported prefix of the Illumina universal adapter. Recovery
requires at least 8 adapter bases and 50 retained read bases. A tentative unique
placement must remain unique at the same strand-aware origin after 8 additional
retained bases are removed.

An accepted recovery preserves the complete input sequence and qualities in
the BAM and represents the adapter tail as a strand-correct soft clip. Its MAPQ
is the lower of the recovery and stability evidence, capped at 20. For a read
that already mapped, only the reported endpoint may change, and only at the
same biological origin; mapping class and MAPQ remain unchanged. Reads without
exact adapter support are unchanged.

!!! note "Single-end confidence"
    Unique origins receive numeric Q10/Q15/Q20/Q30/Q40 from evidence retained
    by the selecting search; tied origins remain MAPQ 0. Output declares
    `caller-compatible-directional-single` or
    `caller-compatible-nondirectional-single`, according to the selected
    library profile. The published single-end truth qualification is scoped to
    directional reads in the documented 5M-R1 simulated corpus.

Single-end output preserves input order. Coordinate-based analysis still
requires coordinate sorting. Add `--non-directional` to run both the OT/OB and
CTOT/CTOB passes and make one global decision. An equal-best cross-pass result
is ambiguous at MAPQ 0; a weaker cross-pass placement contributes to MAPQ
separation and repeat pressure. PBAT remains unsupported.

Shared options including `--sensitive`, `--non-directional`, `--output-contract`,
`--mapped-only`, `--threads`, `--bam-threads`, `--bam-compression-level`, and
`--metrics` apply to the read-1-only layout. Single-end Bismark output takes
the traceback path when auxiliary replay is required; mapped-only output omits
records without an accepted placement. Single-end metrics use the
`bsbit-single-alignment-metrics-v2` schema and include the selected output
policy alongside search work, adapter outcomes, and direct-versus-traceback
record counts. Template-span, total-thread-budget, and paired batching controls
remain paired-only and fail explicitly on single-end input.

## Align paired reads

Use synchronized R1 and R2 for caller-compatible MAPQ and BAM provenance:

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --threads 8 \
  --bam-threads 2 \
  --output-contract minimal \
  --min-template-span 0 \
  --max-template-span 1000
```

On a dedicated host, `--total-threads 14` selects the qualified throughput
split and constrains auxiliary work to its own physical-core pool. A balanced
stride-16 index uses 11 mapping + 3 BGZF workers; a fast stride-8 index uses
10 + 4 because its faster mapping exposes more output pressure. The option
replaces both explicit thread flags:

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --total-threads 14
```

Other budgets reserve approximately one fifth of their cores for output,
bounded to four BGZF workers; fast indexes add one output worker from budgets
of at least 12 when that bound permits. Keep explicit `--threads` and
`--bam-threads` when a scheduler or benchmark requires a particular split.

The default run keeps stdout clean. Add `--metrics` and redirect stdout when
the two-row profiling TSV is useful:

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --threads 8 \
  --metrics \
  > sample.summary.tsv
```

The BAM is written in input order. Inspect it immediately, then follow the
[BAM-preparation guide](prepare-bam.md) before coordinate indexing or calling.

## Choose an alignment mode

Both layouts expose default and sensitive search. For single-end input,
default keeps the existing early-resolution d3/d5 path and `--sensitive`
audits that result against the wider bounded candidate frontier with the Q20
replacement/rescue gate described above. Strong incumbents (Q20 or above and
edit distance at most two) use the sufficient d3 audit boundary; weaker,
distance-three, and unmapped results retain d5 verification. The paired-end
path additionally records one of these stable strategies:

| Mode | Stable strategy | Main behavior |
|---|---|---|
| Default (omit a mode flag) | `balanced-d5-adapter-recovery-read-complete-v2` | Balanced distance-3-to-5 search plus exact Illumina-adapter recovery for otherwise-unmapped directional pairs |
| `--sensitive` | `sensitive-bounded-integrated-read-complete-v1` | Additional complete-frontier, adapter, confidence, and bounded-repeat evidence within the published wall-time envelope |

Use default mode for the normal latency/accuracy balance. Use `--sensitive`
when its higher single-end confidence separation or the documented paired-end
extra recall is worth the additional runtime. Current measurements and their
workload boundary live on the
[performance page](../performance-evidence.md).

Both modes write one primary record per input read by default. An ambiguous
pair may retain one deterministic mapped representative, normally at MAPQ 0;
a narrowly qualified sensitive subset may receive MAPQ 10. Use pair-minimum
MAPQ 20 when a downstream workflow requires pairs classified as unique by this
executable. bsbit MAPQ tiers rank confidence within bsbit and are not
probability-matched to another aligner.

## Select library orientation

Directional libraries are the default for both layouts. Add
`--non-directional` to make one placement decision across all four supported
directional classes. For single-end input, use the command above with this flag;
for paired input, use both mates as shown below. The option does not infer
orientation from assay name and does not change the output-tag contract.

```bash
bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.nondirectional.bam \
  --non-directional \
  --threads 8
```

PBAT is unsupported rather than silently approximated. Check the actual
library protocol before choosing this option.

## Select the output contract

The default `minimal` contract emits `NM` and `XG` and is sufficient for the
in-tree callers. Use `--output-contract bismark` only for a consumer that
requires Bismark-compatible `MD`, `XM`, and `XR` tags in addition to `NM` and
`XG`. The compatibility contract changes tags, not coordinates, ambiguity,
MAPQ, or classification.

See [Alignment BAM](../outputs/alignment-bam.md) for header provenance, record
completeness, tag definitions, and inspection commands.

## Next

- [Prepare the BAM for calling](prepare-bam.md)
- [Review every alignment option](../reference/cli.md#bsbit-align)
- [Review supported workflows](../getting-started/workflow.md#supported-workflows)
- [Troubleshoot alignment](../help/troubleshooting.md)
