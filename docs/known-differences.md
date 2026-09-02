# Known differences and limitations

This page lists differences that affect workflow choice or interpretation.
Historical algorithm candidates are not hidden product modes; the paired-end
CLI exposes only default and `--sensitive`.

## Cross-aligner comparisons

BitMapperBS, HISAT-3N, Bismark, and BISCUIT differ in candidate discovery,
clipping, scoring, ambiguity proof, pair selection, MAPQ, and whether ambiguous
or unmapped reads appear in BAM. Aggregate mapping rate is therefore not an
accuracy oracle, and an integer MAPQ is not necessarily probability-matched
between tools.

The current simulated-truth scorecard, exact tool identities, output-policy
boundary, and contemporaneous runtime comparison are maintained in
[performance evidence](performance-evidence.md). Superseded `strict`,
`semi-global`, and rescue experiment tables are intentionally excluded from
the user documentation.

## Confidence reporting

The caller-compatible paired aligner uses deterministic pair-level score-gap and
repeat evidence, followed in sensitive mode by fixed Q10/Q20/Q30/Q40 operating
tiers. These tiers rank confidence within bsbit. Only their explicitly
documented qualification sets have measured error bounds; they are not a claim
that every record has a universal Phred-calibrated error probability.

Ambiguous pairs may retain one deterministic representative, normally at MAPQ
0. A narrowly bounded ambiguous subset may receive MAPQ 10 but remains
ambiguous. Only native-unique evidence can enter Q20 or higher. A consumer that
requires native-unique paired evidence should therefore use pair-minimum MAPQ
20, not merely `MAPQ > 0`.

Single-end `bsbit align` uses numeric evidence tiers and declares
caller-compatible directional or non-directional provenance. The directional
Q40 tier is qualified within 5 bp on the documented 5M-R1 simulated corpus;
exact-coordinate Q40 does not pass the same gate, and that corpus does not
qualify non-directional single-end.

## Reference and output differences

Reference N bases receive deterministic symbols in the internal search image,
while the exact catalog retains the authoritative N mask.
Candidates cannot cross N or contig barriers, and every accepted placement is
verified against exact reference bases.

`--output-contract bismark` provides compatible `MD`, `XM`, `XR`, and `XG`
evidence for consumers that need those tags. It does not promise byte-identical
Bismark coordinates, CIGAR, flags, MAPQ, ambiguity policy, or clipping. bsbit
does not infer `MM/ML` modification probabilities.

## Unsupported or separately qualified uses

- Directional paired-end WGBS owns the frozen large-corpus performance and
  MAPQ qualification.
- Non-directional paired-end and single-end alignment have four-strand semantic
  and caller-compatibility coverage, not that large-corpus qualification.
- Directional single-end alignment is supported by `bsbit align`, has three
  retained 5M-R1 timing runs, and has exact plus within-5-bp truth evaluation
  on one controlled simulated corpus. That corpus does not qualify
  non-directional reads, other read lengths, assays, references, or hosts.
- Preprocessed RRBS and targeted reads are accepted when chemistry and
  orientation match; bsbit does not trim adapters, model restriction sites,
  interpret targets, or provide assay-specific QC.
- PBAT, CRAM, remote/object-store input, native Windows, macOS, and ARM64 Linux
  are unsupported or unqualified.
- Methylation and SNV callers require study-specific validation. The SNV caller
  is diploid and does not call indels or haplotypes; no clinical-use claim is
  made.

See [sequencing data support](getting-started/workflow.md#sequencing-data-support)
and the [limitations and roadmap](getting-started/workflow.md#limitations-and-roadmap)
for the current scope, and the [scientific contract](scientific-contract.md)
for biological semantics.
