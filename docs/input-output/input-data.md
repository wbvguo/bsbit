# Input data

This page summarizes the inputs required at each stage. See [File
formats](../reference/file-formats.md) for field definitions and examples.

## Inputs by command

| Command | Required inputs | Optional input |
| --- | --- | --- |
| [`bsbit index`](../guides/indexing.md) | Reference FASTA | — |
| [`bsbit align`](../guides/alignment.md) | bsbit reference index and read 1 FASTQ | Read 2 FASTQ for paired-end data |
| [`bsbit call meth`](../guides/methylation-calling.md) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit call snp`](../guides/variant-calling.md) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit call joint`](../guides/variant-calling.md#run-joint-calling) | Prepared BAM and matching reference FASTA | Regions or target BED |
| [`bsbit combine`](../guides/methylation-matrices.md) | Sorted CGmap and/or extended bedMethyl files | — |

Malformed records, unsupported encodings, and incompatible inputs stop the
command with an error.

## Paths and compression

| Input | Accepted form | Sidecars |
| --- | --- | --- |
| FASTA for `index` | Uncompressed, gzip, or BGZF | — |
| FASTA for `call` | Uncompressed or BGZF | Uncompressed: `.fai` recommended; BGZF: `.fai` and `.gzi` required |
| FASTQ | Plain, gzip, or BGZF | — |
| BAM for `call` | Coordinate-sorted BAM | BAI or CSI required |
| Target BED | Plain, gzip, or BGZF | — |
| Methylation calls for `combine` | Plain, gzip, or BGZF | — |

??? note "Local files only"

    Inputs must be regular local files. stdin (`-`), URLs, object-store paths,
    and remote streams are not supported.

??? note "FASTA compression and indexes for `bsbit call`"

    `bsbit call` uses an adjacent `.fai` with an uncompressed FASTA. If the
    index is absent, bsbit issues a warning and scans the FASTA once to build a
    temporary in-memory line-layout index; it does not write a sidecar. A
    BGZF-compressed FASTA requires adjacent `.fai` and `.gzi` indexes. Ordinary
    gzip does not provide random access and cannot be used with `bsbit call`.
    `bsbit index` reads FASTA sequentially and therefore accepts ordinary gzip.

## FASTA reference

- Each record starts with `>` and a nonempty, unique contig name.
- The first whitespace-delimited header token is used as the contig name;
  names are case-sensitive.
- Contig names and lengths must conform to the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
- bsbit accepts at most 1,000,000 contigs and 64,000,000 aggregate contig-name bytes.
- Sequence is case-insensitive and may contain only `A`, `C`, `G`, `T`, and
  `N`.

!!! important "Use the same reference FASTA"

    Use the same reference-genome FASTA for `bsbit index` and every downstream
    `bsbit call` command. If the reference changes, rebuild the bsbit index and
    realign the reads.

## FASTQ reads

bsbit accepts strict, unwrapped four-line FASTQ with 3–192 bases per read.
The first whitespace-delimited header token is used as the read name and must
follow the SAM query-name rules. Sequences may contain only `A`, `C`, `G`, `T`,
and `N`; qualities must use printable Phred+33.

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
checks](../help/troubleshoot.md#paired-fastq-names-or-counts-are-inconsistent)
before a long run.

## Calling BAM

Every `bsbit call` command requires a [prepared BAM](../guides/prepare-bam.md)
created from `bsbit align` output and representing a single biological sample.
If multiple read groups specify `SM`, all nonempty values must match.

Sorting, duplicate marking, and other BAM processing must preserve the standard
per-contig `@SQ M5` reference checksums and mapped-record `XG` tags written by
`bsbit align`. The read layout and library profile recorded in local `@HD bs`
metadata are informational; callers do not branch on them.

## Regions and target BED

Use `--region CONTIG:START-END` to limit calling to a single region. For
multiple regions, use `--regions-bed targets.bed`. `--region` uses 1-based
inclusive coordinates, while BED3+ uses 0-based half-open coordinates. See
[BED target intervals](../reference/file-formats.md#bed-target-intervals) for the accepted
file format.

??? note "Multiple inline regions"

    `--region` accepts `CONTIG:START-END` with 1-based coordinates; both
    `START` and `END` are included. The option can be used more than once in
    the same command and can also be used together with `--regions-bed`. For
    multiple regions, we recommend using a BED file with `--regions-bed`. BED
    uses 0-based, half-open intervals. bsbit merges all supplied intervals
    before calling, so overlaps do not produce duplicate calls.

## Methylation calls for combine

`bsbit combine` accepts one coordinate-sorted CGmap or extended bedMethyl file
per sample. All inputs must use the same reference and compatible coordinates;
the two formats may be mixed. See [Build methylation
matrix](../guides/methylation-matrices.md) for sample naming and merge options.

## Next

- [Build index](../guides/indexing.md)
- [Align reads](../guides/alignment.md)
- [Prepare BAM file](../guides/prepare-bam.md)
- [File formats](../reference/file-formats.md)
- [Troubleshoot rejected input](../help/troubleshoot.md)
