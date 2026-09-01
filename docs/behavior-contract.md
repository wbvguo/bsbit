# Product behavior contract

This page defines behavior that users and downstream tools may rely on. The
[scientific contract](scientific-contract.md) owns bisulfite chemistry and
calling semantics; the [CLI reference](reference/cli.md) owns option spelling
and defaults.

## Command surfaces

`bsbit align` is the only alignment command. Input layout is determined by the
presence of the two explicit read paths:

| Command/layout | Stable role | Input and output |
|---|---|---|
| `bsbit align` with read 1 only | Caller-compatible directional single-end alignment | One FASTQ to BAM in input order; unique origins receive numeric MAPQ |
| `bsbit align` with read 1 and read 2 | Caller-compatible, high-throughput paired-end alignment | Synchronized paired FASTQ to read-complete BAM in input order |

Both layouts accept exactly two search modes: default, selected by omitting a
mode flag, and `--sensitive`. Every other mode spelling is rejected. Paired
profiling summaries record the selected mode and an immutable `strategy_id`;
current identities and measurements live in [performance
evidence](performance-evidence.md).

Directional paired libraries are the default. `--non-directional` performs one
global placement decision across all four supported strand/mate configurations.
PBAT is not silently approximated.

## Input and identity

- FASTA and FASTQ may be plain, gzip, or BGZF local regular files.
- Sequence is normalized case-insensitively to A/C/G/T/N. Malformed records,
  other symbols, mate-name disagreement, and unequal paired EOF are errors.
- Reference contigs and N runs are hard search barriers.
- The exact catalog and internal search data bind to the same semantic
  reference digest. Missing or mismatched bundle data fails before mapping.
- URLs, stdin aliases, devices, and object-store paths are unsupported.

Complete calling requirements, including indexed FASTA and BAM identity, are
defined in [Prepare input data](reference/input-data.md).

## Placement and classification

A single result is `unique` when one best strand-aware biological origin
survives, `ambiguous` when equal-best origins remain, and `unmapped` when no
verified origin survives. A paired result is `unique` when one best concordant
biological placement survives the selected complete policy, `ambiguous` when
equal-best placements remain, and `unmapped` when no accepted pair exists.
Concordance requires the same contig, an admitted library orientation, and a
template span inside the configured inclusive bounds.

Projection is candidate discovery only. Every accepted placement is verified
against the exact four-letter reference under the selected bisulfite strand
relation. Equivalent CIGAR/end-point representations at one strand-aware
5-prime origin count as one placement; equal-score distinct origins remain
ambiguous.

Default mode prioritizes low latency and retries only bounded unresolved work.
For single-end input, `--sensitive` completes the wider six-round, 4,096-hit
candidate frontier before d5 verification and MAPQ. For paired input it also
adds bounded failed-pair, repeat, mate-rescue, and endpoint evidence before
final classification. Resource caps may conservatively retain ambiguity or
leave a read or pair unmapped; they never turn an incomplete frontier into an
unsupported unique claim. Simulator truth, read names, known coordinates, and
peer-aligner output are unavailable to mapping decisions.

After primary classification, an otherwise-unmapped single read with exact
supported 3' Illumina-adapter evidence may be remapped from its retained
prefix. At least 50 bases must remain, and a unique recovery must retain the
same strand-aware origin after an additional 8-base shortening. Accepted BAM
records preserve the complete SEQ and QUAL, use terminal soft clipping, and
cap MAPQ at 20. This fallback never replaces an existing complete-read single
placement.

## BAM output and MAPQ

`bsbit align` writes one primary record per input read by default. Unique pairs
are mapped proper pairs. An ambiguous pair may retain one deterministic mapped
representative, normally at MAPQ 0, while remaining `ambiguous` in the summary.
A pair without a retained placement produces two unmapped primary records.
`--mapped-only` removes only truly unmapped records; it does not mean
`MAPQ > 0` and does not remove mapped MAPQ-0 representatives.

The `minimal` contract emits literal `NM` and conversion-strand `XG`. The
explicit `bismark` contract adds compatible `MD`, `XM`, and `XR` tags without
claiming byte-for-byte Bismark coordinates, classification, or confidence.
Soft-clipped alignments retain the complete original `SEQ` and `QUAL` and use
strand-correct terminal CIGAR operations.

MAPQ is a deterministic within-aligner confidence ranking. Sensitive mode may
assign fixed Q10/Q20/Q30/Q40 operating tiers only to evidence subsets covered
by its policy; an ambiguous representative never enters Q20 or above. These
integers are not universally probability-matched to another aligner. Current
calibration results and their corpus boundary are maintained only on the
[performance page](performance-evidence.md).

Directional single-end `bsbit align` assigns Q10/Q15/Q20/Q30/Q40 from evidence
already retained by the selecting search; MAPQ calculation does not launch a
second search. Tied origins remain MAPQ 0. The BAM declares
`caller-compatible-directional-single` in structured `@PG` provenance and is
accepted by `bsbit call` after the same sorting, indexing, tag, and reference
identity checks as paired output. The current exact and within-5-bp calibration
boundary is reported on the [performance page](performance-evidence.md).

## Determinism and publication

For identical immutable inputs, index bytes, options, worker count, and binary,
classification, record order, and output bytes are deterministic. Paired-end
alignment does not learn an insert-size prior; the explicit template bounds are
the complete span policy.

Every result destination is create-only. Writers stage private bytes and expose
the final target only after successful finalization. Malformed input, reference
identity failure, corrupt dimensions or offsets, arithmetic overflow, resource
failure, worker failure, output collision, and publication failure terminate
without publishing a successful target. A failure never selects a legacy,
high-memory, or otherwise unqualified fallback backend.
