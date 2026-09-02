# Prepare BAM file { #prepare-bam-for-calling }

`bsbit call` requires a coordinate-sorted and indexed BAM. Starting from the
input-order BAM produced by `bsbit align`, prepare it according to the duplicate
policy for the library.

## Decide whether to mark duplicates

PCR or optical duplicates can cause evidence from the same source molecule to
be counted more than once and bias analysis. `samtools markdup` identifies
duplicates from alignment position and orientation, then sets the SAM duplicate
flag so downstream tools can ignore them without removing records. See the
[samtools duplicate-marking
algorithm](https://www.htslib.org/algorithms/duplicate.html) for details.

Coordinate-based duplicate marking is commonly appropriate for randomly
fragmented libraries, but is generally not recommended for RRBS, amplicon, or
other fixed-end libraries because independent molecules may have the same
coordinates by design. UMI libraries require a UMI-aware method. See the
[Bismark deduplication
guidance](https://felixkrueger.github.io/Bismark/usage/deduplication/) for
additional bisulfite-library recommendations.

## Prepare the BAM

If coordinate-based duplicate marking is appropriate for the library, follow
the paired-end or single-end commands below.

**Paired-end data**

```bash
# Name-sort the BAM and add mate information
samtools sort -n -o sample.qname.bam sample.bam
samtools fixmate -m sample.qname.bam sample.fixmate.bam
# Coordinate-sort the BAM, mark duplicates, and index
samtools sort -o sample.sorted.bam sample.fixmate.bam
samtools markdup sample.sorted.bam sample.prep.bam
samtools index sample.prep.bam
```

**Single-end data**

```bash
# Coordinate-sort the BAM, mark duplicates, and index
samtools sort -o sample.sorted.bam sample.bam
samtools markdup sample.sorted.bam sample.prep.bam
samtools index sample.prep.bam
```

`samtools markdup` marks duplicates without removing them. bsbit callers ignore
records marked as duplicates.

??? note "Prepare without duplicate marking"

    When duplicate marking is not part of the selected policy, coordinate-sort
    the alignment BAM and create its index directly:

    ```bash
    samtools sort -o sample.prep.bam sample.bam
    samtools index sample.prep.bam
    ```

## Validate the BAM

Validate the prepared BAM before calling:

```bash
samtools quickcheck -v sample.prep.bam
samtools flagstat sample.prep.bam
```

Sorting and duplicate handling must preserve the bsbit `@PG` header and the
`XG` tags on mapped records. See [Input
data](../reference/input-data.md#calling-bam-and-reference) for the complete
calling requirements and [Alignment BAM](../outputs/index.md#alignment-bam) for BAM
fields and tags.

## Next

- [Call methylation](methylation.md)
- [Call SNVs](variant-calling.md)
- [Input data](../reference/input-data.md#calling-bam-and-reference)
