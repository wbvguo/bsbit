# Architecture

`bsbit` is a Rust workspace with narrow native boundaries. Standard single-end
and caller-compatible high-throughput paired-end alignment share the compact
exact-reference catalog, persisted combined search index, bounded verification
kernels, and bisulfite rules. The layouts differ only where mate pairing,
paired rescue, and calibrated MAPQ require it. The design separates alignment,
indexing, calling, matrix assembly, and publication so each persisted object
can be validated at its owning boundary.

## Data flow

```text
FASTA / FASTA.gz
      |
      v
opaque complete index (.bsbit handle)
      |                                  \
      v                                   v
compact exact-reference catalog     digest-bound combined SA16 image
      |                                   |
      +----------------+------------------+
                       |
             shared candidate search
             and bounded verification
                /              \
               v                v
bsbit align (single)                bsbit align (paired)
one directional FASTQ              directional/non-directional paired FASTQ
      |                                   |
ordered BAM, uncalibrated              ordered BAM
current caller boundary                   |
                           coordinate sort / mark duplicates / BAI or CSI
                                           |
               indexed authoritative FASTA +
                                           v
                                      bsbit call
                 /       |       \
            CGmap   bedMethyl    VCF
                         |
     named bedMethyl samples -> bsbit combine -> matrices
```

The logical index is the semantic authority: contig names, exact bases, N
barriers, and a stable digest. `bsbit index` builds and binds all internal
search data to that digest. Canonical alignment opens the completed bundle
read-only; it never invokes an index builder. Internal search data remains
demand-paged from read-only files. Candidate coordinates from the internal
search image are always checked against the exact catalog
before they become output records.

## Module boundaries

| Module | Responsibility |
|---|---|
| `bsbit-core` | Stable values and invariants only: bases, normalized sequences, semantic reference identity, bisulfite chemistry, coordinates, and structural CIGAR |
| `bsbit-io` | Format-neutral file lifecycle: absolute-path validation, create-new staging, identity checks, synchronization, and atomic publication |
| `bsbit-hts` | Biological file formats and transport: FASTA, FASTQ, BED3, the shared SAM/BAM alignment model, SAM text, BAM/HTSlib, BGZF text, and shared bedMethyl records |
| `bsbit-index` | Reference ownership plus index construction and storage: the packed reference catalog and projected FM/rank/locate search image |
| `bsbit-align` | The complete alignment domain: edit distance, CIGAR replay, scalar/SIMD kernels, seeds, candidates, extension, pairing, paired-end modes and phases, rescue/search policy, and MAPQ |
| `bsbit-call` | Alignment-evidence reconstruction, fragment overlap collapse, regional aggregation, methylation/SNV likelihoods, and call-specific CGmap/VCF rendering |
| `bsbit-combine` | Parallel preflight and memory-efficient merge of named bedMethyl inputs into matrices |
| `bsbit-cli` | Product commands, cross-crate composition, runtime policy, and the single executable |
| `crates/bsbit-hts/htslib-shim` | Small project-owned C ABI over pinned HTSlib |

Dependencies point from product composition toward reusable mechanisms and are
checked by `tests/tools/test_crate_boundaries.py`:

```text
bsbit-cli     -> bsbit-call, bsbit-combine, bsbit-align,
                 bsbit-index, bsbit-hts, bsbit-io, bsbit-core
bsbit-call    -> bsbit-core, bsbit-hts, bsbit-io
bsbit-combine -> bsbit-hts, bsbit-io
bsbit-align   -> bsbit-index, bsbit-core
bsbit-index   -> bsbit-io, bsbit-core
bsbit-hts     -> bsbit-io, bsbit-core
bsbit-io      -> std/platform
bsbit-core    -> std, SHA-256 mechanism
```

The arrows are dependency ceilings, not a reason to add an otherwise unused
dependency. In particular, `bsbit-io` has no biological-format dependency;
`bsbit-hts` has no alignment or index dependency; and `bsbit-align` has no HTS
dependency. Converting an alignment result and reference dictionary into a
shared SAM/BAM alignment record therefore happens at the `bsbit-cli` composition
boundary. Tests may use lower-level crates as development dependencies to
construct fixtures without changing the dependency graph.

Within `bsbit-cli`, `command/align.rs` owns the one user-visible alignment
entry point. It selects single-end input when only `--read1` is present and
paired-end input when both `--read1` and `--read2` are present. The CLI owns
FASTQ/BAM orchestration and record composition; both alignment algorithms live
in `bsbit-align`.

Callers and matrix assembly use `bsbit-io` directly only for format-neutral
path validation. Biological decoding, encoding, and BGZF transport continue
to flow through `bsbit-hts`.

For high-throughput output, `crates/bsbit-hts/src/alignment_record.rs` also owns the
compact borrowed alignment-record contract and its worker-local batch. The
application can retain those records without first constructing the owned
general record; both the BAM field encoder and the SAM encoder consume that
same contract, with the SAM writer resolving reference ordinals through its
header dictionary.
The BAM-specific binary CIGAR encoding and HTSlib calls remain in `bam.rs`.

The practical boundary test is based on inputs and outputs, not a word in a
symbol name. Comparing an already selected reference slice with a query and
returning an edit distance, endpoint, or CIGAR belongs to
`bsbit_align::verification`.
Code that accepts a reference index, seed/candidate state, mate constraints, or
MAPQ policy belongs to `bsbit_align::search` or the flat alignment-domain
modules. FM rank/locate operations and reference ownership belong to
`bsbit-index`; persistent index protocols belong to `bsbit_index::storage`;
their constructors belong to `bsbit_index::build`. Stable DNA, coordinate, and
semantic reference-identity values belong to `bsbit-core`.

Directory depth is deliberately limited. A crate keeps stable top-level
concepts as files and introduces one directory level only for a cohesive
family such as `build`, `storage`, `verification`, `search`, `paired_end`,
`meth`, `snp`, `evidence`, `region`, or `command`.
Files are split on ownership and change boundaries, not on a line-count quota;
a tightly coupled implementation may remain large rather than expose private
state merely to make the tree look symmetrical.

Tracked `src/` contains only selected product behavior. New algorithms, policy
switches, mapping modes, and ablations begin in dated `agent/worktree/` attempts. A
successful attempt enters a crate under a clear domain name; rejected or
superseded evidence remains with that attempt instead of becoming a dormant
Cargo feature or a copied `experiments/` source tree.

Safe crates forbid unsafe Rust. Architecture intrinsics, raw syscalls, and
native FFI are isolated in the modules that document their safety invariants.
Fixed third-party source is under `external/` and is authenticated before use.

## Data and persistence invariants

The `ReferenceSemanticDigest` defined in `bsbit-core` is the format-neutral
semantic authority that binds every derived product. The opaque index handle
is a packed reference catalog containing validated contig names, exact bases,
aggregate dimensions, and per-contig integrity checks; its private sibling
search image is bound to the same digest. Indexing and calling depend
independently on the core digest value; the caller has no dependency on
`bsbit-index`. Contig identifiers,
contig-local offsets, global offsets, FM rows, query offsets, and SAM positions
remain distinct checked domains rather than interchangeable integers.

The paired search surface combines a demand-paged BWT, exact 16-mer lookup,
occurrence checkpoints, and SA16 rank samples with a separately validated
packed reference catalog. Search projections may collapse letters, but every
accepted placement is verified against the exact bases, contig boundaries, and
N mask owned by that catalog. A projected domain larger than `u32` requires a
new 64-bit format; it is never truncated.

Mapping workers reuse bounded per-worker state for projected queries, FM
intervals, located rows, candidate windows, verification, traceback, and pair
selection. Calling workers keep bounded indexed-reference windows and only the
mate observations that may overlap. Matrix merging retains one current row per
input plus bounded merge channels, so missing cells remain distinct from zero
without accumulating every genomic site in memory.

Persisted product formats use fixed little-endian fields, explicit
magic/version values, checked lengths, and authenticated content. Writers stage
immutable bytes and publish the visible target or metadata last. Readers verify
identity before exposing mapped slices and reject unknown versions, invalid
ordering, overflow, truncation, trailing data, and digest mismatches. The index
is intentionally opaque at the command boundary; its user-visible artifact and
publication contract is documented in the [outputs overview](outputs/index.md).

## Calling flow

`bsbit-call` is the internal calling-library boundary. Its `meth`, `snp`, and
`joint` modules expose independent typed options and `call` entry points. The
user-facing `bsbit call` command is parsed in `bsbit-cli` and delegates the
scientific work to this crate. The internal crate name does not create a second
user-facing executable.

Mode-specific aggregation, likelihood, rendering, and orchestration live under
`meth`, `snp`, and `joint`. BAM fragment reconstruction is owned once by
`evidence`; region selection and planning are owned by `region`. Input
preflight, indexed-reference context, bounded region workers, and create-only
publication remain explicit shared files rather than a catch-all calling
directory.

The caller projects observed bases from BAM CIGAR and SEQ onto the required
indexed reference FASTA; `MD` is ignored. Required `XG`, FLAG, base-quality, and
mapping-quality fields supply the remaining evidence. R1/R2 observations from one
fragment are matched by QNAME, read group, and reciprocal coordinates, then
collapsed at overlapping positions by active-filter eligibility, canonical
base status, and combined base/mapping error, with deterministic R1
tie-breaking. Only potentially overlapping records are
cached. Up to 64 genomic lanes are classified through two-bit base planes and
accumulated in bit-sliced counters. Adaptive indexed coordinate regions provide
bounded parallelism; joint calling shares the first-pass evidence between
methylation and SNV candidate detection.
The public caller is single-sample: the BAM header may declare multiple read
groups only when their `SM` values agree. Optional user targets are normalized
from direct intervals and BED3+ into one dictionary-ordered union before those
internal work regions are created.
`bsbit-hts` owns the format boundary: indexed BAM access, decoded SEQ/CIGAR/aux
fields, BGZF encoding, and format-aware finalization. Generic staging,
identity checks, synchronization, create-only publication, and rollback are
delegated to `bsbit-io`.
Completed region results are rendered directly to staging output in ordinal
order through bounded channels and a sliding reorder window, so the caller
never accumulates a whole-genome site vector. Worker panics are caught and
reported through the same structured error channel. SNV region sizing accounts
for worst-case retained candidates and calls as well as bit-sliced counters;
exact likelihood batches have a separate per-worker planning budget. At caller
startup, indexed FAI/GZI-backed FASTA is streamed in bounded chunks to verify
the BAM's semantic reference digest. Region workers then fetch only the spans
needed for authoritative methylation context and SNV reference alleles. This
applies to `meth`, `snp`, and `joint`.

## Index construction

`bsbit index` invokes the bounded internal search-index builder after creating
the exact reference catalog. It projects the validated catalog into the three-symbol search domain,
constructs exact libsais32 blocks, merges them in memory, and publishes
digest-bound metadata last. The image contains a dense exact 16-mer lookup,
BWT, Occ64/Occ65536 checkpoints, and SA16 samples. Unknown bases receive
deterministic projection symbols, while the catalog's N mask remains
authoritative.

Publication is create-only, bundle-atomic, and fail-closed. Structural checks
are mandatory; the public command has one fixed bounded construction path.
`bsbit align` contains no construction path and fails when required internal
data is absent.

## Alignment flow

`bsbit align` opens the opaque index read-only and chooses the input layout
explicitly from the supplied read paths. With only `--read1`, it runs the
deterministic directional single-end alignment and preserves FASTQ order. With
both `--read1` and `--read2`, it runs the paired-end path described below. Both
layouts publish BAM through the same create-only output contract.

Single-end default mode may classify a verified initial frontier immediately
and continues only unresolved work. Single-end sensitive mode first preserves
that exact default result as its incumbent, then replays intervals admitted by
the shared 4,096-hit sensitive bound, completes the six-round adaptive seed
schedule, and verifies the accumulated candidates through distance 5. A
different-origin replacement or new rescue is accepted only when the completed
frontier is unique at MAPQ 20 or above. Lower-confidence conflicts retain the
default representative, classify it as ambiguous, and emit MAPQ 0; an
uncertified rescue remains unmapped. The single-end path does not enter
mate-rescue, template-geometry, or paired adapter stages.

One decode stream reads single-end records; paired input uses two synchronized
decode streams. Both feed bounded batches to mapping workers with reusable
workspaces. Completed results pass through a bounded ordered queue to record
composition and BAM output, which preserve input order and publish the final
file only after successful finalization.

Strict whole-read mapping is always the first pass. With explicit windowed mate
rescue, a verified uniquely seeded mate can constrain four exact proof-block
queries for its missing partner to the configured local pair window. The
resulting sparse frontier uses the ordinary flexible verifier; it does not
materialize the whole-reference word-parallel filter planes. In default mode,
a global proof interval above the locate cap may take a second complete path:
one unique best low-distance anchor constrains a `d + 1`-block rolling scan to
the same mate window. All starts still pass the flexible verifier and complete
tie selection; no global interval is silently truncated.

Sensitive mode adds a distinct bounded failed-pair stage after the default
path. It ranks each mate's `d + 1` exact-block intervals by rarity, locates at
most 512 cumulative global rows, and performs the pair-geometry join before
the d5 verifier. If only one mate is informative, its verified placements
drive the same bounded-window completion used by rescue. Because this global
frontier can be partial, the workspace propagates a completeness bit into pair
selection: a recovered result from a partial frontier is downgraded to
ambiguous and has its second-best claim removed. Its serialized representative
starts at MAPQ 0 and may enter only the bounded empirical Q10 tier; it cannot
enter Q20/Q30/Q40. Both public modes emit deterministic ambiguous
representatives and preserve pairs without a representative as unmapped
primary records.

The same ranked-block completion is also applied to a nominally unique result
when window rescue contributed or either mate located at least 256 rows. This
is a bounded second-best search: a newly found lower-score pair replaces the
original, and a distinct equal-score origin makes the result ambiguous. Search
and confidence thresholds are deliberately separate. After completion, a
sensitive unique result that still depends on window rescue or has at
least 384 located rows is capped at MAPQ 19.

After this bounded search, sensitive mode sends only structural ambiguous
pairs (a retained candidate clips bases or changes reference/query span) to a
worker-owned, narrow-band affine scorer. Conversion-aware BWA score units
rerank that residual set; already unique pairs never pay this dynamic-program
cost. MAPQ uses the BWA paired score-gap and near-best-count form, while every
remaining ambiguity is serialized deterministically with baseline MAPQ 0.
The aligner does not apply an insert-size prior; explicit template-span bounds
are the only span policy.

The final sensitive reporting layer applies fixed integer confidence
certificates after the BWA-style baseline and repeat/clipping caps. A bounded
ambiguous subset may enter Q10 while retaining native ambiguous class;
native-unique complete-frontier evidence supplies the additional Q30 and Q40
tiers. These comparisons are allocation-free and use only counters and
geometry already produced by alignment. Truth is used only by the ignored
development qualification, never by alignment.

Only after that proof and affine work does sensitive mode consider semi-global
endpoint search. The gate admits a complete, non-unmapped frontier below MAPQ
20. The candidate must keep the same strand-aware pair origin and increase the
effective evidence MAPQ; otherwise the original result is restored. This
prevents clipping from discovering the first coordinate or laundering a
partial search frontier into a high-confidence call.

In the fixed adapter-trimmed phase, complete-read unmapped pairs with qualified
3' adapter evidence enter a compact second batch. A
window-rescue ambiguity whose nominal baseline had no pair geometry is also
eligible, so rescue cannot suppress an otherwise valid clipped recovery.
Tentative unique recoveries use one additional shortened batch to require a
stable strand-aware 5' origin. Original read bytes remain owned by the input
batch and are serialized in full; the direct BAM layer adds terminal soft-clip
CIGAR operations without creating a separately trimmed FASTQ representation.

Sensitive mode's qualified endpoint phase stays inside each worker's existing
strict candidate frontier and uses the targeted gate above. It computes
conversion-aware mismatch and barrier prefixes once per nominal placement and
chooses bounded 5'/3' retained endpoints with a linear terminal scan. The
strict and clipped placement sets are globally re-sorted before the ordinary
spatial pair join. Equal-scoring CIGAR representations sharing the same
strand-aware mate 5' origins collapse to one biological placement; distinct
origins remain ties. A retained distance-zero tentative unique receives a
complete FM repeat check, with oversized repeat anchors failing closed to
ambiguous. Ungapped clipped output uses a direct slab-to-BAM path and does not
run traceback merely to reconstruct `S/M/S`.

`--threads` controls mapping workers; the operating system schedules those
workers normally. Decode, ordered record construction, BAM writing, queues,
and batch sizes have explicit bounded resource contracts. The aligner does not
perform per-read development audits or read repository-local data.

Directional mode runs the ordinary paired read-conversion configuration.
`--non-directional` runs both the ordinary and complementary configurations as
complete searches and merges their best, second-best, and near-best evidence
before classification and MAPQ. A cross-configuration tie remains ambiguous.
This library selection is independent of the public default and sensitive
search modes.

## Persistence and failure model

Every product format has magic/version fields, checked lengths, fixed
endianness, and digests or structural validation. Paths are create-only.
Truncation, checksum mismatch, unsupported format, arithmetic overflow,
resource-limit violation, or publication-identity change is an error; no
backend silently falls back to a less qualified implementation.
