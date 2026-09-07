# Alignment generalization protocol v1

This protocol separates policy development from evidence used to support
scientific claims. Scenarios not yet executed are requirements for the final
manuscript validation matrix, not claims about current results.

## Frozen decision rule

- Alignment policy: `bounded-structural-alignment-v1`
- Alignment policy source: `crates/bsbit-align/src/alignment_policy.rs`
- MAPQ policy: `structural-origin-evidence-v1`
- MAPQ policy source: `crates/bsbit-align/src/mapq_policy.rs`
- Validation begins only after both policy identifiers and all numeric fields are
  frozen.
- Any threshold or decision-rule change increments the identifier and
  invalidates prior qualification results for the changed policy.
- Refactors may retain the identifier only when exhaustive unit/regression
  tests demonstrate identical classifications, coordinates, and MAPQ values.

No truth coordinate, simulator metadata, benchmark tool output, or dataset
label may enter runtime alignment or MAPQ decisions. Reference ordinals may
address index storage, but they must never act as an otherwise-equal reporting
preference or MAPQ feature.

## Development diagnostics

Candidate-frontier traces are excluded from normal builds. Enable the explicit
`development-trace` feature only for offline audits:

```console
cargo run -p bsbit-cli --release --features development-trace --example candidate_trace -- ...
cargo run -p bsbit-cli --release --features development-trace --example paired_stage_trace -- ...
```

The feature may retain bounded candidate and verified-mate details and is not a
benchmark or production setting. It must remain disabled for all reported time,
memory, throughput, mapping-rate, and accuracy measurements.

## Dataset separation

Use three disjoint sets:

1. **Development:** may be inspected while designing features and thresholds.
2. **Pre-registered validation:** scenario parameters and seeds are committed
   before results are opened; no policy changes are allowed afterward.
3. **External confirmation:** independently generated or real spike-in data
   not used during development.

Report development and validation results separately. Never pool them into one
headline accuracy value.

## Required factorial coverage

Every available scenario is run as SE and PE, directional and non-directional,
and default and sensitive.

| Axis | Required held-out levels |
| --- | --- |
| Reference architecture | Mammalian primary assembly; alternate human assembly; compact microbial genome; repeat-enriched synthetic panel |
| Read length | 50, 100, and 150 bases, plus the production default |
| Substitution error | Low, nominal, and high profiles generated with independent seeds |
| Indel burden | None/low and elevated profiles |
| Bisulfite conversion | 95%, 98%, and 99.5% conversion efficiency |
| Methylation composition | Low, mixed, and high methylation; CG and non-CG contexts where applicable |
| PE insert distribution | Short, nominal, and long means with narrow and broad dispersion |
| Repeat stratum | Unique sequence, low-copy repeats, and high-copy/segmental duplication sequence |
| Coverage sample | At least five independent read-generation seeds per scenario family |

Reference families and truth-generation seeds used during policy development
must not be reused as the sole evidence for the validation set.

## Comparators and execution controls

- Compare bsbit default and sensitive with the same external aligner versions
  used in the performance matrix.
- Record complete command lines, binary hashes, index hashes, input hashes,
  host CPU, thread allocation, wall time, CPU time, peak RSS, and output policy.
- Use the same allowed library direction, template-span knowledge, thread
  budget, and read-output accounting for every tool.
- Distinguish output records from reads or pairs with valid mapped coordinates.
- Evaluate coordinates before filtering by MAPQ, then evaluate cumulative Q10,
  Q20, Q30, and Q40 subsets.

## Primary endpoints

For every scenario and mode report:

- output-read fraction and mapped fraction;
- exact-coordinate and within-five-base correctness;
- precision, recall, and F1 at Q10, Q20, Q30, and Q40;
- empirical error rate among records passing each threshold;
- one-sided 95% Clopper--Pearson upper error bound;
- wall time, CPU time, peak RSS, and reads/s or pairs/s.

A Q threshold passes calibration only when the one-sided 95% upper error bound
is at most `10^(-Q/10)`. Empty threshold sets are reported as `N/A`, never as
precision 1.

## Pre-registered acceptance rules

The frozen policy is considered to generalize only when:

1. Every claimed Q threshold passes its calibration requirement in the pooled
   validation set and in each declared reference family.
2. Sensitive mode has the highest or statistically indistinguishable F1 among
   calibrated tools in a majority of scenarios, with Q10--Q30 designated as
   primary endpoints and Q40 as a high-confidence endpoint.
3. No individual scenario shows an unexplained concentrated failure in a read
   length, conversion rate, repeat stratum, or insert-size stratum.
4. Default-mode accuracy remains within the pre-registered non-inferiority
   margin while retaining its intended performance advantage.
5. Results are reproduced from clean indexes and immutable inputs by a second
   execution manifest.

If a validation result motivates a policy change, move that scenario into the
development set, increment the policy identifier, and register a new untouched
validation set before making a new generalization claim.
