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

For a library using coordinate-based duplicate marking, name-sort the BAM, add
mate information, coordinate-sort it, mark duplicates, and create the BAM
index:

```bash
samtools sort -n -o sample.qname.bam sample.bam
samtools fixmate -m sample.qname.bam sample.fixmate.bam
samtools sort -o sample.sorted.bam sample.fixmate.bam
samtools markdup sample.sorted.bam sample.prep.bam
samtools index sample.prep.bam
```

`samtools markdup` marks duplicate records rather than removing them; bsbit
callers ignore records with the duplicate flag.

??? note "Prepare without duplicate marking"

    When duplicate marking is not part of the selected policy, coordinate-sort
    the alignment BAM and create its index directly:

    ```bash
    samtools sort -o sample.prep.bam sample.bam
    samtools index sample.prep.bam
    ```

## Validate the BAM

Check the final BAM before calling:

```bash
samtools quickcheck -v sample.prep.bam
samtools flagstat sample.prep.bam
```

The prepared BAM must retain the bsbit `@PG` header entry and `XG` tags. Calling
also requires the same reference FASTA used for alignment. See the [calling
input requirements](../reference/input-data.md#calling-bam-and-reference) for
the complete validation contract and [Alignment BAM](../outputs/alignment-bam.md)
for BAM fields, tags, and provenance.

## Next

- [Call methylation](methylation.md)
- [Call SNVs](variant-calling.md)
- [Review the calling input requirements](../reference/input-data.md#calling-bam-and-reference)
