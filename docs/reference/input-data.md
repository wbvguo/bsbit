# Prepare input data

bsbit reads local FASTA and FASTQ files using strict validation. Malformed
records and unsupported encodings are reported as errors.

`bsbit align` accepts one FASTQ or synchronized paired FASTQ for
caller-compatible directional or explicitly non-directional alignment. Neither
layout requires genome-wide coverage:
preprocessed WGBS, RRBS, and targeted bisulfite read sets share these input
contracts. That compatibility does not add RRBS-specific trimming,
restriction-site handling, target-panel interpretation, or assay-specific QC
to bsbit.

## Compression and paths

FASTA and FASTQ may be plain, gzip, or BGZF. Compression is decoded directly;
do not pre-decompress FASTQ as a performance workaround.

Reference-index construction and alignment stream sequentially through FASTA
and FASTQ files. Calling requires random access: FASTA needs an adjacent `.fai`
(and `.gzi` for BGZF), while BAM needs BAI or CSI. Ordinary gzip-compressed
FASTA is not supported for calling.

Inputs must be regular local paths. stdin (`-`), URLs, object-store paths, and
remote streaming are unsupported.

## Calling BAM and indexed reference

Every `bsbit call meth`, `bsbit call snp`, and `bsbit call joint` run requires:

- a coordinate-sorted BAM with an adjacent BAI or CSI;
- mapped primary observations with an available MAPQ value; MAPQ 255 records
  are excluded rather than interpreted as high confidence;
- one unique structured bsbit `@PG` header line (program/run metadata, not a
  per-read tag) declaring a caller-compatible alignment mode and exact reference
  semantic digest, plus a string `XG:Z:CT|GA` tag on every mapped primary
  record;
- an authoritative FASTA with an adjacent `.fai`, plus `.gzi` when that FASTA
  is BGZF-compressed; and
- one biological sample per BAM. Multiple read groups are accepted only when
  all nonempty `SM` fields agree.

The FASTA must be the same assembly used for alignment. The caller fetches it
in BAM dictionary order, normalizes case, recomputes the semantic digest, and
compares it with BAM provenance. A same-name, same-length FASTA with different
bases is rejected. The caller projects `SEQ` through `CIGAR` and ignores `MD`.

Current single- and paired-end `bsbit align` outputs declare caller-compatible
provenance and satisfy their documented MAPQ contracts. MAPQ 255 remains an
unavailable score if encountered in external or older BAMs and is never
interpreted as high confidence.

Follow [BAM preparation](../guides/prepare-bam.md#prepare-bam-for-calling) before
running a caller. Choose `meth`, `snp`, or `joint` in the [command-line
reference](cli.md#choose-a-command).

## FASTA reference

- A record begins with `>` and a nonempty name.
- The name is the first whitespace-delimited token in the header.
- Names are case-sensitive and must be unique in the reference.
- Sequence is normalized to uppercase.
- The accepted alphabet is `A`, `C`, `G`, `T`, and `N` (case-insensitive).
- Other IUPAC ambiguity codes are rejected rather than silently collapsed.

Use exactly the same FASTA assembly and contig naming convention for every
artifact in one alignment workflow.

## FASTQ reads

FASTQ is strict four-line input:

```text
@read-name
ACGTT...
+
IIIII...
```

- The record name is the first whitespace-delimited header token.
- Sequence must be nonempty and use `A`, `C`, `G`, `T`, or `N`.
- Sequence and quality must each occupy one line and have equal lengths.
- Quality bytes must be printable Phred+33 characters.
- If the `+` line repeats a header suffix, it must agree with the record header.

## Paired-read synchronization

The paired command takes mates through `bsbit align --read1 R1 --read2 R2`.
R1 and R2 must end together and remain synchronized by ordinal and name.
Accepted name forms are:

- identical source names, such as `instrument:run:read`; or
- the matching pair `instrument:run:read/1` and
  `instrument:run:read/2`.

The shared stem becomes the BAM query name. Reordered files, missing records,
or inconsistent names stop the run.

When provenance is uncertain, use the [paired-FASTQ troubleshooting
checks](../help/troubleshooting.md#paired-fastq-names-or-counts-are-inconsistent)
before a long run. Counts alone do not prove name synchronization; alignment
validates every name and pair.

## Library orientation

`bsbit align` supports directional single-end or paired-end reads by default.
With `--non-directional`, either layout makes one placement decision across all
four supported bisulfite directions. It does not silently approximate PBAT. The
exact conversion/orientation relations are defined in the [scientific
contract](../scientific-contract.md).

Do not infer library orientation from coverage design. WGBS, RRBS, and targeted
bisulfite experiments may each produce libraries with different orientation or
preprocessing requirements; compatibility is determined by the actual library
protocol, not by how much of the reference was assayed.

## Template span

Standard paired `bsbit align` uses inclusive `--min-template-span` and
`--max-template-span` bounds, defaulting to 0 and 1000. These paired-only flags
are rejected for standard single-end input. Pairs may overlap or one mate may
contain the other when the selected bounds permit it.

## Next

- [Build index](../guides/indexing.md)
- [Align reads](../guides/alignment.md)
- [Troubleshoot rejected input](../help/troubleshooting.md)
