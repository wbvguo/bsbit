# Input data

This page summarizes the inputs required at each stage. See [File
formats](file-formats.md) for field definitions and examples.

## Inputs by command

| Command | Required inputs | Optional input |
|---|---|---|
| [`bsbit index`](../guides/indexing.md) | Reference FASTA | — |
| [`bsbit align`](../guides/alignment.md) | bsbit reference index and read 1 FASTQ | Read 2 FASTQ for paired-end data |
| [`bsbit call meth`](../guides/methylation.md) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit call snp`](../guides/variant-calling.md) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit call joint`](../guides/variant-calling.md#run-joint-calling) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit combine`](../guides/methylation-matrices.md) | Sorted CGmap and/or extended bedMethyl files | — |

Malformed records, unsupported encodings, and incompatible inputs stop the
command with an error.

## Paths and compression

| Input | Accepted form | Sidecars |
|---|---|---|
| FASTA for `index` | Plain or BGZF | None |
| FASTA for `call` | Plain or BGZF | Plain: `.fai` optional; BGZF: `.fai` and `.gzi` required |
| FASTQ | Plain, gzip, or BGZF | None |
| BAM for `call` | BAM | Adjacent `.bai` required |
| Target BED | Plain, gzip, or BGZF | None |
| Methylation calls for `combine` | Plain, gzip, or BGZF | None |

Inputs must be regular local files. stdin (`-`), URLs, object-store paths, and
remote streams are not supported.

## FASTA reference

- Each record starts with `>` and a nonempty, unique contig name.
- The first whitespace-delimited header token is used as the contig name;
  names are case-sensitive.
- Sequence is case-insensitive and may contain only `A`, `C`, `G`, `T`, and
  `N`.

!!! important "Use the same reference FASTA"

    Use the same reference-genome FASTA for `bsbit index` and every downstream
    `bsbit call` command. If the reference changes, rebuild the bsbit index and
    realign the reads.

## FASTQ reads

bsbit accepts strict, unwrapped four-line FASTQ. The sequence must be nonempty,
contain only `A`, `C`, `G`, `T`, or `N`, and have the same length as the
printable Phred+33 quality string. If the `+` line repeats the header, the two
must agree.

Preprocessed WGBS, RRBS, and targeted bisulfite reads follow the same rules.
bsbit does not perform assay-specific trimming or read QC.

## Paired-read synchronization

Supply mates with `bsbit align -1 R1 -2 R2`. The files must contain the same
number of records in the same order. Accepted name forms are:

- identical names; or
- matching `/1` and `/2` suffixes.

The shared name becomes the BAM query name. Reordered reads, missing mates, or
inconsistent names stop the run.

When synchronization is uncertain, use the [paired-FASTQ troubleshooting
checks](../help/troubleshooting.md#paired-fastq-names-or-counts-are-inconsistent)
before a long run.

## Library orientation

Directional alignment is used by default. Add `--non-directional` for a
non-directional single-end or paired-end library; bsbit then makes one
placement decision across all four supported bisulfite directions.

Choose the setting from the library protocol, not from whether the experiment
is WGBS, RRBS, or targeted. See the [alignment guide](../guides/alignment.md)
for usage and the [scientific contract](../development/scientific-contract.md) for the
conversion model.

## Calling BAM and reference

Every `bsbit call` command requires:

- a coordinate-sorted BAM with an adjacent `.bai`;
- one caller-compatible bsbit `@PG` header and an `XG:Z:CT|GA` tag on every
  mapped primary record;
- numeric MAPQ values; records with MAPQ 255 are excluded; and
- one biological sample per BAM. If multiple read groups define `SM`, their
  nonempty values must agree.

The reference must be the same FASTA used for alignment. The caller verifies
its normalized sequence digest against the BAM provenance, so matching contig
names and lengths alone is not sufficient.

Follow [Prepare BAM file](../guides/prepare-bam.md#prepare-bam-for-calling)
before calling.

## Regions and target BED

Calling can be limited with either `--region CONTIG:START-END` or
`--regions-file targets.bed`, but not both. Inline regions use 1-based inclusive
coordinates; BED uses 0-based half-open coordinates. See [BED target
intervals](file-formats.md#bed-target-intervals) for the accepted file format.

## Methylation calls for combine

`bsbit combine` accepts one coordinate-sorted CGmap or extended bedMethyl file
per sample. All inputs must use the same reference and compatible coordinates;
the two formats may be mixed. See [Build methylation
matrix](../guides/methylation-matrices.md) for sample naming and merge options.

## Next

- [Build index](../guides/indexing.md)
- [Align reads](../guides/alignment.md)
- [Prepare BAM file](../guides/prepare-bam.md)
- [File formats](file-formats.md)
- [Troubleshoot rejected input](../help/troubleshooting.md)
