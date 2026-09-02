# Scientific behavior contract

This page defines the biological and coordinate semantics that alignment and
calling must preserve. The [product behavior contract](behavior-contract.md)
owns CLI-visible classification, output completeness, MAPQ, determinism, and
publication behavior.

## Evidence and scope

The bisulfite model follows the chemistry described by
[Frommer et al.](https://doi.org/10.1073/pnas.89.5.1827), including the fact
that conventional bisulfite sequencing cannot distinguish 5mC from 5hmC
without additional chemistry, as discussed by
[Huang et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2811190/). Strand names
follow the library vocabulary documented by
[Bismark](https://felixkrueger.github.io/Bismark/usage/library-types/).

Directional paired WGBS is the large-corpus qualified alignment surface.
Explicit non-directional paired and single-end alignment have four-strand
semantic and compatibility coverage. Directional single-end alignment has one
controlled large-corpus same-binary default/sensitive performance and truth
comparison, but neither single-end profile inherits the replicated paired
performance/MAPQ qualification, and the directional single-end comparison does
not qualify non-directional reads. PBAT, EM-seq, TAPS,
oxBS-seq, TAB-seq, hairpin bisulfite sequencing, long-read modification
calling, and protocol-specific RRBS processing require separate contracts and
are not silently approximated.

Adapter/primer preprocessing, duplicate policy, conversion-rate QC, and
study-level biological interpretation remain pipeline responsibilities.

## Sequence and coordinates

- A/C/G/T are canonical. Input is normalized to uppercase for computation
  without changing the retained read record.
- `N` is unknown, never a wildcard or a zero-cost match. Other IUPAC symbols
  are rejected.
- Contigs and N runs are hard candidate barriers.
- Internal intervals are 0-based half-open. SAM `POS` conversion occurs only
  at serialization.
- Contig IDs, global and local offsets, query positions, FM rows, and template
  spans remain checked domains. Overflow is an error.

## Bisulfite strand relation

After a read is oriented left-to-right in forward-reference coordinates, the
supported strand identities are:

| Strand | SAM orientation | Cytosine evidence | Zero-cost conversion |
|---|---|---|---|
| OT | forward | reference `C` | reference `C`, query `T` |
| CTOT | reverse | reference `C` | reference `C`, query `T` |
| OB | reverse | reference `G` | reference `G`, query `A` |
| CTOB | forward | reference `G` | reference `G`, query `A` |

For canonical bases, the exact relations are:

```text
top:    match(r, q) = (r == q) OR (r == C AND q == T)
bottom: match(r, q) = (r == q) OR (r == G AND q == A)
```

The relation is asymmetric: reference T/query C and reference A/query G are
ordinary substitutions. Retained and converted observations both have zero
conversion-aware substitution cost at the relevant cytosine. Mapping must not
prefer a retained base merely because it is literally equal to the reference.

A retained C or G is evidence of a protected or unconverted cytosine, not proof
of 5mC. Incomplete conversion, 5hmC, sequencing error, variation, and mapping
error remain possible explanations.

## Alignment invariants

C-to-T and G-to-A projections discover candidates; they cannot determine final
distance, CIGAR, ambiguity, or methylation because projection collapses
biologically distinct symbols. Every accepted placement is verified against
the exact four-letter reference under the selected strand relation.

Insertion, deletion, and non-bisulfite substitution cost one in the complete-
read edit model. Terminal clipping is admitted only by the documented bounded
recovery policy. It never rewrites the original read: BAM retains complete
`SEQ` and `QUAL`, and the CIGAR identifies the clipped portion.

Pair selection compares biological strand-aware 5-prime origins, not merely
serialized coordinates. Equivalent CIGAR/end-point representations of one
origin collapse to one placement; equal-best distinct origins remain
ambiguous. A bounded heuristic may decline an unresolved frontier or retain
ambiguity, but it cannot certify unique from incomplete evidence. Runtime
alignment has no simulator truth, expected coordinate, read-name exception, or
cross-aligner oracle.

Scalar, SSE4.2, AVX2, and any available AVX-512 verification kernels must
produce identical distance, endpoint, and tie sets. Traceback uses the same
frozen scoring and tie policy. Scheduling, worker count, and batch partitioning
must not change classification or record order.

## Library profiles

For single-end reads, directional mode searches OT and OB. Non-directional mode
also searches CTOT and CTOB, merges evidence before classification, and retains
an equal-best cross-pass result as ambiguous.

Directional mode admits the first two template classes. Non-directional mode
also admits the complementary read orders:

| Template class | R1 | R2 | Inward order |
|---|---|---|---|
| original top | OT / forward | CTOT / reverse | R1 is not right of R2 |
| original bottom | OB / reverse | CTOB / forward | R2 is not right of R1 |
| complementary top | CTOT / reverse | OT / forward | R2 is not right of R1 |
| complementary bottom | CTOB / forward | OB / reverse | R1 is not right of R2 |

Concordant mates must share a contig, satisfy a class admitted by the selected
profile, and have an outer span inside the configured inclusive bounds.
Equality, overlap, and containment are valid when those rules hold.
Non-directional evidence is merged before classification; an equal-best tie
between configurations remains ambiguous.

## SAM/BAM semantics

Serialization follows the
[SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[standard tags specification](https://samtools.github.io/hts-specs/SAMtags.pdf).

- Reverse alignments serialize reverse-complemented `SEQ`, reversed `QUAL`,
  and orientation-consistent CIGAR and mate fields.
- Complete-read output uses `M/I/D`; bounded terminal recovery may add `S`.
  A bisulfite-compatible conversion is not represented as literal `=`.
- `NM` counts literal differences from the forward reference, not the internal
  conversion-aware edit distance.
- `XG:Z:CT|GA` records the genome-conversion strand. Together with FLAG
  orientation it preserves OT/OB/CTOT/CTOB identity.
- The opt-in Bismark contract adds `MD` and `XM/XR`. `XM` follows stored BAM
  sequence orientation; insertions and soft clips are `.`, deletions consume
  no query byte, and unknown context is reported rather than guessed.
- `MM/ML` modification probabilities are not inferred from bisulfite evidence.

Tag compatibility does not change placement, ambiguity, MAPQ, flags, or
coordinates. See the [behavior contract](behavior-contract.md#bam-output-and-mapq)
for record completeness and confidence semantics.

## Downstream calling semantics

Every `meth`, `snp`, and `joint` run requires a coordinate-sorted BAI/CSI-
indexed BAM with the canonical bsbit `@PG` record and the authoritative indexed
FASTA. Contig names and lengths must agree, and the normalized reference
semantic digest must equal the digest stored by alignment. The caller projects
`SEQ` through `CIGAR` onto FASTA bases and ignores `MD`; a different
same-dictionary FASTA is rejected by the digest check.

Every mapped primary record must contain `XG:Z:CT|GA`. Unmapped, secondary,
supplementary, QC-failed, duplicate-marked, and MAPQ-255 records do not
contribute. One BAM represents one biological sample; multiple read groups are
valid only when their nonempty `SM` fields agree.

Overlapping mates are matched by query name, read group, and reciprocal
coordinates. One fragment contributes at most one observation per genomic
site. Filter eligibility, canonical/present/known-quality status, combined base
and mapping error, and finally R1 provide deterministic tie-breaking shared by
methylation and SNV evidence.

At forward-reference C, retained C is methylated evidence and T is
unmethylated evidence on the top conversion strand. At forward-reference G,
retained G and converted A supply the corresponding bottom-strand evidence.
CG, CHG, and CHH context comes from the authoritative FASTA, including flanks
outside a read or requested region. The operational label “methylated” is not
proof of 5mC, is not corrected for conversion efficiency, and cannot
distinguish 5mC from 5hmC. Opposite members of a CpG dyad remain separate
genomic sites.

SNV discovery excludes directly conversion-confounded top C-to-T and bottom
G-to-A changes, then applies the documented depth and ALT thresholds. Exact
calling evaluates unordered diploid genotypes with base and mapping quality,
strand-specific conversion rates, a reference-divergence prior, and adaptive
integration over unknown methylation for C/G genotypes. `AQ` and `GQ` remain
distinct. The caller emits SNVs only; it is not an indel, haplotype, or
clinical caller.

`joint` shares the overlap-collapsed first evidence pass. `bsbit combine`
performs no new biological inference: it preserves methylated/total counts,
represents absent or filtered cells as missing rather than zero, and filters by
the configured valid-sample proportion. Exact schemas live in the
[methylation](guides/methylation.md),
[variant](guides/variant-calling.md), and
[matrix](guides/methylation-matrices.md) guides.

## Required validation

The maintained tests cover the canonical-base relation and N behavior;
coordinate and reverse-complement round trips; scalar/SIMD equivalence;
FM/rank/locate agreement with naive references; directional, paired, repeat,
indel, overlap, and contig-boundary fixtures; corruption rejection;
deterministic BAM fields; clipping and ambiguity behavior; authoritative-FASTA
calling with overlap collapse; regional equivalence; and deterministic matrix
merging with missing-cell semantics.

Support boundaries are summarized under
[sequencing data support](getting-started/workflow.md#sequencing-data-support). Measured accuracy and
MAPQ evidence live in [performance evidence](performance-evidence.md), and
remaining differences are listed in [known limitations](known-differences.md).
