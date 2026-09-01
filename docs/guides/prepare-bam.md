# Prepare a BAM for calling { #prepare-bam-for-calling }

`bsbit align` writes caller-compatible single- or paired-end records in input order. Inspect
the new BAM before transforming it, but do not try to index it until it is
coordinate-sorted:

```bash
samtools quickcheck -v sample.bam
samtools flagstat sample.bam
```

The paired-end recipe below produces the coordinate-ready
`sample.analysis.bam` used throughout the calling guides.

## Prepare paired-end output

Group mates, populate the mate tags used for duplicate decisions,
coordinate-sort, mark duplicates, and create the BAM index:

```bash
samtools sort -n -o sample.name.bam sample.bam
samtools fixmate -m sample.name.bam sample.fixmate.bam
samtools sort -o sample.position.bam sample.fixmate.bam
samtools markdup sample.position.bam sample.analysis.bam
samtools index sample.analysis.bam
samtools quickcheck -v sample.analysis.bam
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
insufficient. The FASTA needs `.fai` plus `.gzi` when it is BGZF-compressed.

The shared [calling input contract](../reference/input-data.md#calling-bam-and-indexed-reference)
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

- [Call methylation](../outputs/methylation.md)
- [Call SNVs](../outputs/variant-calling.md)
- [Review the calling input contract](../reference/input-data.md#calling-bam-and-indexed-reference)
