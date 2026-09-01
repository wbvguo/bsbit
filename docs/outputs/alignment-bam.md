# Alignment BAM

`bsbit align` writes a SAM/BAM 1.6 stream in input order. The BAM preserves the
alignment decision, reference identity, read sequence and qualities, and the
bisulfite strand information required by bsbit callers. It is not
coordinate-sorted when first written.

For the full artifact sequence, see the [outputs overview](index.md). To use a
paired-end BAM for calling, follow [Prepare a BAM for calling](../guides/prepare-bam.md).

## Header and provenance

The header contains ordinary `@HD` and `@SQ` records plus exactly one canonical
bsbit program record. In SAM text it has this shape:

```text { .no-copy }
@PG ID:bsbit PN:bsbit VN:VERSION DS:reference-semantic-sha256=DIGEST;alignment-mode=MODE
```

The digest identifies the normalized reference sequence rather than only its
path, contig names, or lengths. `MODE` records one of the supported alignment
contracts:

| Mode | Meaning for downstream calling |
|---|---|
| `caller-compatible-directional-single` | Directional single-end BAM with numeric MAPQ, accepted after preparation and reference checks |
| `caller-compatible-nondirectional-single` | Four-strand single-end BAM with numeric MAPQ, accepted after the same preparation and checks |
| `caller-compatible-directional-paired` | Directional paired-end BAM accepted after preparation and reference checks |
| `caller-compatible-nondirectional-paired` | Four-strand paired-end BAM accepted after the same preparation and checks |

Sorting and duplicate-marking tools must preserve this `@PG` record. Adding
their own program records is normal.

## Record completeness and ordering

By default, paired alignment writes one primary record for every input read:

- unique pairs are mapped proper pairs;
- an ambiguous pair may retain a deterministic mapped representative, normally
  at MAPQ 0, while remaining classified as ambiguous;
- a pair without a retained placement produces two unmapped primary records;
  and
- `--mapped-only` removes only truly unmapped records, not mapped MAPQ-0
  representatives.

Single-end input also writes one primary record per input read. Unique origins
receive numeric evidence tiers; tied mapped representatives use MAPQ 0 and
unmapped reads remain explicit.

The initial header declares unsorted order and records follow FASTQ input
order. Run `samtools sort` before coordinate indexing; do not create a BAI or
CSI directly for the original alignment BAM.

## Standard fields

Mapped records use ordinary SAM fields and conventions:

| Field | bsbit contract |
|---|---|
| `FLAG` | Encodes mapping, mate, proper-pair, and reverse-complement state |
| `RNAME`, `POS`, `CIGAR` | Describe placement against the original forward reference |
| `RNEXT`, `PNEXT`, `TLEN` | Describe the paired placement when a mapped mate is present |
| `MAPQ` | Within-aligner confidence tier; not a cross-aligner probability scale |
| `SEQ`, `QUAL` | Preserve the complete read; reverse mappings use SAM orientation |

Bisulfite-compatible conversions remain literal differences in `SEQ`; they
are not rewritten to match the reference. CIGAR uses `M/I/D`, with `S` only
for admitted bounded terminal recovery.

## Auxiliary tag contracts

Choose the tag set with `--output-contract`:

| Contract | Mapped-record tags | Use it when |
|---|---|---|
| `minimal` (default) | `NM`, `XG` | Using bsbit callers or a consumer that needs standard alignment plus conversion-strand identity |
| `bismark` | `NM`, `XG`, `MD`, `XM`, `XR` | A downstream consumer explicitly requires Bismark-compatible methylation and conversion tags |

The tags mean:

| Tag | Meaning |
|---|---|
| `NM:i` | Literal differences from the forward reference, not the internal conversion-aware edit distance |
| `XG:Z:CT|GA` | Genome-conversion strand; required on every mapped primary record used by bsbit callers |
| `MD:Z` | Canonical reference-difference string, emitted only by the `bismark` contract |
| `XM:Z` | Bismark-compatible per-query methylation/context string in stored BAM sequence orientation |
| `XR:Z:CT|GA` | Bismark-compatible read-conversion identity |

Selecting `bismark` changes only the auxiliary tag set. It does not change
placement, coordinates, ambiguity, MAPQ, flags, or classification, and it does
not claim byte-for-byte equivalence with Bismark output. bsbit callers project
`SEQ` through `CIGAR` against the authoritative FASTA and ignore `MD`.

## Inspect alignment output

```bash
samtools quickcheck -v sample.bam
samtools view -H sample.bam
samtools flagstat sample.bam
samtools view sample.bam | head
```

No output from `samtools quickcheck` means its structural checks passed. Use
`samtools view -H` to verify the bsbit provenance before any downstream tool
rewrites the header.

## Caller boundary

The current caller requires caller-compatible single- or paired-end provenance,
mapped primary records with available MAPQ and `XG`, a coordinate-sorted
BAI/CSI-indexed BAM, and the exact indexed FASTA used by alignment. Both
directional and non-directional single-end modes satisfy that provenance
boundary.

See the [calling input contract](../reference/input-data.md#calling-bam-and-indexed-reference)
for the complete validation rules and the [scientific contract](../scientific-contract.md#sambam-semantics)
for strand, coordinate, and methylation semantics.
