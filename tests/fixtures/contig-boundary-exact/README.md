# Contig-boundary exact-search fixture

This small regression proves that exact search cannot create a biological hit
by concatenating adjacent contigs without a barrier.

- `reference.fa` contains two 256-base, C-free contigs.
- `cross-boundary.fastq` contains one 75-base read.
- The read equals `ctgA:[219,256) || ctgB:[0,38)` in zero-based half-open
  coordinates.
- It occurs zero times in either contig and exactly once only in the invalid
  separator-free concatenation.
- The C-free construction makes C-to-T projection irrelevant.

| File | SHA-256 |
|---|---|
| `reference.fa` | `6d5e9b0b0b5350941e4ca4c1be3b79170975f80e04cd86cb9726254ee2488412` |
| `cross-boundary.fastq` | `202460d05661f3cfbf8c940c73aaf8a9c37a5689d5b9b39f4b88ed375197fd0d` |

`crates/bsbit-index/tests/reference.rs` embeds both files and asserts the
contig-local barrier behavior. The fixture is source data for that automated
test; it is not a benchmark input or a compatibility oracle.
