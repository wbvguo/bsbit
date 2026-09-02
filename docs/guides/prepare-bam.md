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

Duplicate marking is a study-design decision: identical starts can be expected
in amplicon and other fixed-end designs. Select a protocol-appropriate policy
rather than removing reads automatically. bsbit callers ignore records
carrying the SAM duplicate flag.

The final BAM must retain bsbit's structured `@PG` provenance and `XG` tags.
`@PG` is one standard SAM/BAM header line describing the program and run
contract, not a tag repeated on every read. It records the exact normalized
reference fingerprint and whether the alignment mode is caller-compatible.
Calling recomputes the normalized semantic digest of the FASTA and requires it
to equal the BAM digest; matching contig names and lengths alone is
insufficient. A plain FASTA uses an existing adjacent FAI when available and is
otherwise scanned once into an in-memory position table. A BGZF-compressed
FASTA requires both `.fai` and `.gzi`; ordinary gzip FASTA is unsupported.

The shared [calling input contract](../reference/input-data.md#calling-bam-and-reference)
lists all checks in one place. The [alignment BAM reference](../outputs/alignment-bam.md)
defines record completeness, MAPQ, provenance, and auxiliary tags.

## Prepare single-end output

Single-end output does not need mate repair. Coordinate-sort it, apply the
study-appropriate duplicate policy when required, and create the BAM index:

```bash
samtools sort -o sample.analysis.bam sample.bam
samtools index sample.analysis.bam
samtools quickcheck -v sample.analysis.bam
```

If the study requires duplicate marking, insert a validated single-end
duplicate workflow before indexing; do not apply paired-only mate-repair
assumptions to single-end records.

Current output declares `caller-compatible-directional-single` or
`caller-compatible-nondirectional-single` and carries numeric MAPQ, so the
caller accepts it after the same reference, `XG`, sort, and index checks. Older
`standard-directional-single` BAMs used unavailable MAPQ 255 and must be
realigned rather than relabeled.

## Next

- [Call methylation](methylation.md)
- [Call SNVs](variant-calling.md)
- [Review the calling input requirements](../reference/input-data.md#calling-bam-and-reference)
