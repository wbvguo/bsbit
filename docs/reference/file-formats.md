# File formats

This page defines the formats, coordinates, and fields used by bsbit. See
[Input data](input-data.md) for accepted inputs and [Output
files](../outputs/index.md) for workflow artifacts.

## Quick reference

| Artifact | Typical name | Used by | Format and transport |
|---|---|---|---|
| [Reference sequence](#fasta-reference) | `reference.fa` | `index`, `call` | FASTA; plain or BGZF |
| [Sequencing reads](#fastq-reads) | `sample_R1.fastq.gz` | `align` | Strict four-line FASTQ; plain, gzip, or BGZF |
| [Target intervals](#bed-target-intervals) | `targets.bed.gz` | `call` | BED3 or BED3+; plain, gzip, or BGZF |
| [bsbit alignment index](#bsbit-alignment-index) | `reference.bsbit` | `index`, `align` | Opaque index handle |
| [BAM alignments](#bam-alignments-and-index) | `sample.bam` | `align`, `call` | SAM/BAM 1.6 in BGZF-compressed BAM |
| [BAM index](#bam-alignments-and-index) | `sample.bam.bai` | `call` | BAI random-access sidecar |
| [CGmap calls](#cgmap-methylation-calls) | `sample.cgmap.gz` | `call meth`, `call joint`, `combine` | Eight-column tab-delimited CGmap |
| [Extended bedMethyl](#extended-bedmethyl-calls) | `sample.bed.gz` | `call meth`, `call joint`, `combine` | Eighteen-column bedMethyl |
| [Variant calls](#vcf-variant-calls) | `sample.vcf.gz` | `call snp`, `call joint` | VCF 4.3 |
| [Methylation matrix](#methylation-matrices) | `cohort.level.bed.gz` or `cohort.count.bed.gz` | `combine` | BED6 plus sample columns |
| [Alignment metrics](#alignment-metrics-tsv) | `alignment.summary.tsv` | `align --metrics` | Two-row profiling TSV written to stdout |

## Coordinate systems

Check the coordinate convention before moving values between formats:

| Format or field | Start/position convention | End convention |
|---|---|---|
| SAM/BAM `POS` | 1-based | Derived from `CIGAR` |
| `--region CONTIG:START-END` | 1-based | Inclusive |
| CGmap column 3 | 1-based | One genomic base |
| VCF `POS` | 1-based | One base for bsbit SNVs |
| BED3, extended bedMethyl, and matrix BED columns 2–3 | 0-based | Half-open; end is excluded |

For example, CGmap position `101` and the BED interval `100 101` identify the
same base.

## FASTA reference

FASTA contains one or more named reference sequences:

```text
>chr1
ACGTCGATCGATCG
>chr2 optional description
NNACCGTT
```

- A record starts with `>` and a nonempty header. The first
  whitespace-delimited token is the contig name, so the second record above is
  named `chr2`.
- Contig names are case-sensitive and must be unique.
- Sequence is case-insensitive on input and is normalized to uppercase.
- bsbit accepts only `A`, `C`, `G`, `T`, and `N`; other IUPAC ambiguity codes
  are rejected.

See [FASTA input requirements](input-data.md#fasta-reference) for compression,
sidecars, and reference-consistency rules.

## FASTQ reads

bsbit accepts strict, unwrapped four-line FASTQ:

```text
@read0001/1 optional description
ACGTTGCA
+
IIIIIIII
```

| Line | Content |
|---:|---|
| 1 | `@` followed by the read name; the first whitespace-delimited token is used |
| 2 | Nonempty sequence containing only `A`, `C`, `G`, `T`, or `N` |
| 3 | `+`, optionally followed by a header suffix that agrees with line 1 |
| 4 | Printable Phred+33 quality characters, exactly one per sequence base |

Sequence and quality wrapping is not supported. See [paired-read
synchronization](input-data.md#paired-read-synchronization) for mate-name and
ordering requirements.

## BED target intervals

`--regions-file` reads the first three tab-separated BED columns. Additional
columns are allowed and ignored:

```text
# capture panel
chr1	100	150	exon-1
chr2	500	575	exon-2
```

| Column | Meaning |
|---:|---|
| 1 | Contig name, which must exist in the BAM dictionary |
| 2 | 0-based start |
| 3 | Exclusive end |
| 4+ | Optional uninterpreted annotations |

Blank lines, `#` comments, and UCSC `track` or `browser` directives are
ignored. Intervals must be nonempty and within the referenced contig when used
for calling. A regions file may be plain, gzip, or BGZF.

## bsbit alignment index

`bsbit index` produces an opaque index with no editable text representation or
public field schema. Treat it as one complete bundle: do not modify its
components, and rebuild it when the reference changes.

## BAM alignments and index

`bsbit align` writes SAM/BAM 1.6 records as BAM. The initial BAM follows FASTQ
input order and therefore must be coordinate-sorted before its `.bai` index is
created. This is a simplified SAM-text view of the binary file:

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1000
@PG	ID:bsbit	PN:bsbit	VN:VERSION	DS:reference-semantic-sha256=DIGEST;alignment-mode=caller-compatible-directional-paired
read0001	99	chr1	101	40	8M	=	181	88	ACGTTGCA	IIIIIIII	NM:i:1	XG:Z:CT
```

The first eleven record fields are standard SAM fields:

| Field | Meaning in bsbit output |
|---|---|
| `QNAME`, `FLAG` | Read name and SAM bit flags |
| `RNAME`, `POS`, `CIGAR` | Placement against the forward reference |
| `MAPQ` | Within-aligner confidence tier |
| `RNEXT`, `PNEXT`, `TLEN` | Mate placement and template length |
| `SEQ`, `QUAL` | Complete stored read sequence and base qualities |

Mapped records always carry `NM:i` and `XG:Z:CT|GA`. With
`--output-contract bismark`, they also carry `MD:Z`, `XM:Z`, and `XR:Z:CT|GA`.
The structured `@PG` line binds the BAM to the exact reference and alignment
mode and must survive sorting and duplicate marking.

`.bai` is a binary random-access sidecar, not an alignment file. A caller
requires a coordinate-sorted BAM and its adjacent `.bai`. See [Alignment
BAM](../outputs/index.md#alignment-bam) for output behavior and validation, and
[Prepare BAM file](../guides/prepare-bam.md) for sorting and indexing.

## CGmap methylation calls

CGmap has no header. Each tab-delimited row describes one cytosine-strand site:

```text
chr1	C	101	CG	CG	0.750000	3	4
chr1	G	102	CG	CG	0.500000	2	4
```

| Column | Meaning |
|---:|---|
| 1 | Contig |
| 2 | Forward-reference base: `C` for the `+` cytosine strand or `G` for the `-` strand |
| 3 | 1-based position |
| 4 | Context: `CG`, `CHG`, or `CHH` |
| 5 | Cytosine-strand dinucleotide |
| 6 | Methylated fraction to six decimal places, or `na` when total coverage is zero |
| 7 | Methylated observation count |
| 8 | Total valid methylated plus unmethylated coverage |

Rows are sorted by BAM dictionary, coordinate, and strand. Column 7 cannot
exceed column 8, and matrix construction derives the level from the counts
rather than trusting a rounded value in column 6.

## Extended bedMethyl calls

The extended bedMethyl output uses exactly 18 tab-separated columns and
0-based, half-open one-base intervals:

```text
chr1	100	101	m,CG,0	4	+	100	101	255,0,0	4	75.00	3	1	0	0	0	0	0
```

| Column | Meaning | Example |
|---:|---|---|
| 1 | Contig | `chr1` |
| 2 | 0-based start | `100` |
| 3 | Exclusive end; start + 1 | `101` |
| 4 | Modification and context: `m,CG,0`, `m,CHG,0`, or `m,CHH,0` | `m,CG,0` |
| 5 | Valid coverage | `4` |
| 6 | Cytosine strand: `+` or `-` | `+` |
| 7 | Display start; equal to column 2 | `100` |
| 8 | Display end; equal to column 3 | `101` |
| 9 | Display color | `255,0,0` |
| 10 | Valid coverage; equal to column 5 | `4` |
| 11 | Percent methylated to two decimal places | `75.00` |
| 12 | Methylated observations | `3` |
| 13 | Unmethylated observations | `1` |
| 14 | Other-modification observations | `0` |
| 15 | Deletion observations | `0` |
| 16 | Failed observations | `0` |
| 17 | Different-base observations | `0` |
| 18 | No-call observations | `0` |

Columns 5 and 10 must both equal columns 12 + 13. Deletions and different-base
observations are reported separately and do not enter valid coverage. `bsbit
call meth` reports zero for other modification, failed, and no-call counts.

CGmap is compact and compatible with CGmap-oriented tooling. Extended
bedMethyl carries explicit strand, interval, and evidence categories and is
directly compatible with the BED coordinate model. Both can be mixed as input
to `bsbit combine`. See [Call methylation](../guides/methylation.md) for the
evidence rules shared by both representations.

## VCF variant calls

`bsbit call snp` and the SNV half of `bsbit call joint` write VCF 4.3. A file
contains `##` metadata, one `#CHROM` header, and one row per called SNV. This
excerpt omits most metadata declarations but shows the complete data columns:

```text
##fileformat=VCFv4.3
##source=bsbit
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
chr1	101	.	A	G	50.25	PASS	DP=8;AF=0.500000;BSI=ONE;BS8=4,0,2,0,2,0,0,0	GT:GQ:AQ:DP:AD:IAD:PL	0/1:42:37:8:6,2:4,2:50,0,50
```

| Fixed field | Meaning |
|---|---|
| `CHROM`, `POS` | Contig and 1-based SNV position |
| `ID` | `.`; bsbit does not assign variant identifiers |
| `REF`, `ALT` | Reference base and one or more alternate bases |
| `QUAL` | Confidence that the site is not homozygous reference |
| `FILTER` | `PASS`, `LowAD`, `LowGQ`, `LowAQ`, or a semicolon-separated combination |
| `INFO` | Site-level depth, allele fraction, and bisulfite-strand evidence |
| `FORMAT`, sample | Names and values of the one biological sample |

Important bsbit fields are:

| Field | Location | Meaning |
|---|---|---|
| `DP` | INFO/FORMAT | Quality-filtered fragment depth |
| `AF` | INFO | Conditional expected ALT fraction |
| `BSI` | INFO | Whether discrimination uses `BOTH` strands or `ONE` unaffected strand |
| `BS8` | INFO | Top-strand A,C,G,T counts followed by bottom-strand A,C,G,T counts |
| `GT` | FORMAT | Unphased maximum-likelihood genotype dosage |
| `GQ` | FORMAT | Complete diploid-genotype confidence |
| `AQ` | FORMAT | Confidence that each selected ALT is present |
| `AD`, `IAD` | FORMAT | Raw and bisulfite-informative allele depths |
| `PL` | FORMAT | Prior-free normalized genotype likelihoods |

bsbit writes SNVs only, not indels. See [Call SNVs](../guides/variant-calling.md)
for filtering and statistical interpretation.

## Methylation matrices

`bsbit combine` writes four metadata lines, a BED6-plus-samples header, and
one row per retained site. A level matrix has one fraction per sample:

```text
##bsbit_matrix_format=level
##bsbit_min_count=1
##bsbit_min_prop=0.000000000
##bsbit_cg_only=false
#chrom	start	end	modification	score	strand	tumor	normal
chr1	100	101	m,CG,0	0	+	0.750000	0.500000
chr1	101	102	m,CG,0	0	-	.	0.250000
```

A count matrix has paired methylated and total-count columns:

```text
##bsbit_matrix_format=count
##bsbit_min_count=1
##bsbit_min_prop=0.000000000
##bsbit_cg_only=false
#chrom	start	end	modification	score	strand	tumor_meth_count	tumor_total_count	normal_meth_count	normal_total_count
chr1	100	101	m,CG,0	0	+	3	4	2	4
chr1	101	102	m,CG,0	0	-	.	.	1	4
```

The first six columns use the extended bedMethyl coordinate, modification, and
strand model. `.` means an absent or below-threshold sample cell; it never
means numeric zero. With `-m both`, the level and count schemas are
written to separate `.level` and `.count` files. See
[Build methylation matrix](../guides/methylation-matrices.md) for filtering,
naming, and mixed-input normalization.

## Alignment metrics TSV

`bsbit align --metrics` writes a self-describing two-row TSV to stdout. Redirect
it separately from the BAM:

```bash
bsbit align \
  -x reference.bsbit \
  -1 reads.fastq.gz \
  -o sample.bam \
  --metrics \
  > sample.alignment.summary.tsv
```

The first row is the column header. The second row starts with a
layout-specific versioned schema identifier:

- `bsbit-single-alignment-metrics-v2` for single-end alignment; or
- `bsbit-alignment-metrics-v2` for paired-end alignment.

Columns cover:

- read or pair counts and BAM record counts;
- mapping, BAM, library, search, and output settings;
- decode, mapping, queue, compression, and output-finalization timings;
- soft-clip fallback and mate-rescue counts; and
- MAPQ, strategy, and read-output policy identifiers.

Metrics are profiling diagnostics, not alignment records or caller input.
