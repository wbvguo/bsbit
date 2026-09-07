# Alignment algorithm and mapping confidence

This page defines the scientific alignment contract implemented by `bsbit
align`. It describes observable decision rules rather than architecture-specific
SIMD kernels or worker scheduling.

The candidate-search and endpoint policy identifier is
`bounded-structural-alignment-v1`; its frozen numeric table lives in
`crates/bsbit-align/src/alignment_policy.rs`. MAPQ is versioned independently.

## Implementation ownership

| Concern | Owning module |
| --- | --- |
| Versioned search, adapter, endpoint, and proof bounds | `alignment_policy.rs` |
| SE/PE search orchestration and structural evidence collection | `single_end/mapper.rs`, `paired_end/mapper.rs` |
| Origin grouping and cross-pass evidence merge | Layout-specific `mapper/merge.rs` |
| Reference-order-independent reporting lottery | `reporting_tie_break.rs` and PE `mapper/reporting.rs` |
| All confidence calculations, caps, certificates, and MAPQ-trigger decisions | Layout-specific `mapq.rs` |
| Versioned numeric MAPQ table only | `mapq_policy.rs` |
| Bounded offline candidate diagnostics | Feature-gated `mapper/trace.rs` |

Production mapper and merge modules do not read `MAPQ_POLICY` directly. This
keeps evidence generation separate from confidence interpretation and makes the
frozen rule table auditable without tracing search control flow.

## Terms

- A **placement** is one verified query/reference endpoint representation.
- A **biological origin** groups placements that describe the same molecule
  origin despite equivalent indel shifts, clipping endpoints, or bisulfite
  projections.
- A **frontier** is the bounded set of candidates and alternatives examined by
  one search policy.
- A **certificate** is truth-blind evidence collected from the reference and
  read that permits a confidence declaration. It never uses a simulated truth
  coordinate.

## Alignment pipeline

For each read or read pair, alignment follows the same conceptual stages:

```text
for each conversion pass allowed by the library profile:
    project the query into the combined bisulfite alphabet
    run adaptive exact offset-seed searches against the combined FM-index
    merge repeated seed evidence by candidate coordinate
    verify candidates with conversion-aware bounded edit distance
    group equivalent endpoints into biological origins

merge conversion-pass results by the same origin-level objective
if paired-end:
    enforce inward orientation and the configured template-span interval
    attempt bounded mate rescue when the selected policy permits it
    rerank the residual qualified frontier with affine endpoint scores
assign MAPQ from retained origin evidence
select a reproducible representative for an unresolved MAPQ-0 tie
```

Candidate generation and verification rank genomic evidence. The reporting
tie-break is applied only after classification and MAPQ are fixed, so changing
the tie-break seed cannot turn an ambiguous result into a unique result.

## Library profiles

| Profile | Conversion passes | Final decision |
| --- | --- | --- |
| Directional | Original | One origin-level decision over OT and OB labels |
| Non-directional | Original and complementary | One joint decision after merging OT, OB, CTOT, and CTOB origin evidence |

The complementary paired-end pass swaps the input mates for the canonical
directional executor, relabels its molecular strands, and restores input mate
order before the two pass results are merged.

## Search modes

Both modes use the same index, conversion-aware edit objective, maximum edit
distance, origin definition, output contract, alignment-policy version, and
MAPQ policy. They differ in how much evidence is collected before a decision
is finalized.

| Component | Default | Sensitive |
| --- | --- | --- |
| Common verification | Distance-three first pass with incremental distance-five fallback | Same common path retained as an incumbent |
| Adaptive seed hit budget | 1,000 per seed | 4,096 per seed |
| SE offset-seed rounds | Up to 6 | Up to 10 |
| PE offset-seed rounds | Up to 6 | Up to 6 plus ranked disjoint-block completion |
| Provisional unique audit | Narrow corroboration of eligible default SE results | Complete bounded audit of SE results and eligible PE frontiers |
| Mate rescue | Disabled except adapter-stability handling | Bounded window rescue from verified anchors |
| Endpoint completion | Adapter-supported fallback | Bounded semi-global endpoint search with affine reranking |
| High-confidence evidence | May be limited by an incomplete frontier | Requires the applicable completed-frontier certificate |

“Complete” always means complete with respect to the documented bounded
frontier, not exhaustive dynamic programming against every reference base.

The run-level endpoint policy is shared by SE and PE. In `auto` mode, exact
3′ adapter evidence can trigger a clipped stability remap in either layout,
and sensitive PE can additionally complete endpoints within an already
discovered candidate locus. `adapter` retains only the exact adapter path;
`none` disables both paths. Endpoint completion cannot introduce an
undiscovered five-prime origin, and its total query clipping is bounded by
`--max-soft-clip`.

## Origin selection

The primary objective is conversion-aware whole-read edit evidence. Paired-end
selection additionally requires compatible molecular strands, inward
orientation, and an outer template span inside the configured inclusive
interval. When multiple endpoint representations belong to one biological
origin, only that origin's best endpoint contributes to runner-up evidence.
This prevents indel shifts or alternative clipping representations at one
locus from being counted as independent genomic competitors.

Sensitive affine scoring is a secondary objective on a bounded residual
frontier. It does not replace candidate discovery and cannot create a genomic
origin that was absent from the verified frontier.

## Structural MAPQ

The policy identifier is `structural-origin-evidence-v1`. MAPQ is an
evidence-derived structural confidence score, not a fitted posterior
probability.

For single-end alignment, the baseline is

```text
raw MAPQ = min(60, 10 × whole-read edit separation)
```

where the separation is the best/runner-up distance difference. If no runner-up
was observed, the next unverified distance outside the completed verification
boundary supplies the conservative alternative.

For paired-end alignment, the baseline uses the best/runner-up origin score
gap, with ten MAPQ units per score unit and a logarithmic multiplicity penalty
for additional near-best origins.

The final decision is deliberately ordered:

1. **Raw separation:** compute confidence from observed origin separation.
2. **Adverse caps:** incomplete search, repeat pressure, high edit burden,
   rescue provenance, instability, clipping, or an observed near-best origin
   may only lower confidence.
3. **Positive certificates:** a small named set of completed-frontier,
   independent-seed certificates may support a declared confidence floor or
   boundary.

### Certificate registry

| Certificate | Required observed evidence | Use |
| --- | --- | --- |
| Completed multi-seed coordinate | At least three independent offset seeds rediscover one coordinate in a completed frontier | Coordinate corroboration |
| Long-read singleton corroboration | A direct singleton on a read of at least 128 bases is independently rediscovered, with edit distance at most two | Coordinate corroboration |
| Long-read boundary corroboration | Independent long-read seeds support a winner at the completed verification boundary with no runner-up | Coordinate corroboration |
| Low-edit local-locus corroboration | Two offsets support one already-collapsed origin with edit distance at most one | Local endpoint equivalence |
| Audited short-singleton coordinate | A bounded 16--46 base suffix becomes a singleton and a bounded high-confidence audit completes | High-confidence origin evidence |
| Audited strong multi-seed coordinate | Multiple offsets support the same low-edit winner and a bounded high-confidence audit completes | High-confidence origin evidence |
| Completed pair alternative margin | The pair frontier is extended beyond the winner sufficiently to expose close alternative origins | High-confidence pair evidence |

Certificate names describe reference/read evidence and deliberately do not
contain a target Q threshold. The numerical relationship between a certificate
and an output tier is versioned separately in the policy table.

Default mode performs the bounded confidence audit only for provisional unique,
low-edit, high-confidence SE results.  A completed default audit may therefore
support MAPQ 40 or greater without claiming completion of the wider sensitive
candidate frontier.  Sensitive mode completes that wider frontier for every
read.  An unaudited default result remains subject to the uncertified
high-confidence cap.

## Ambiguous output and reproducibility

An unresolved equal-best set remains ambiguous with MAPQ 0. By default bsbit
still writes one primary record for every input read. When at least one verified
origin exists, its representative is selected by a seed-controlled hash of the
read identity, contig name, molecular strand, and five-prime coordinate. The
lottery is reproducible and independent of reference contig order.

## Calibration claim

The implementation alone does not assert that an integer score is exactly
`-10 log10 P(error)`. Q10, Q20, Q30, and Q40 operating points must be reported
with held-out empirical precision, recall, F1, error rate, and a one-sided 95%
binomial upper confidence bound. A threshold is called calibrated only when
that upper bound does not exceed the error probability implied by the declared
Q value.
