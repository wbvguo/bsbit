<img src="docs/img/bsbit.png" alt="bsbit logo" width="150" align="right">

# bsbit

bsbit is an ultrafast, memory-efficient toolkit for bisulfite-sequencing
analysis. It provides reference indexing, read alignment, methylation and
single-nucleotide variant calling, and cohort-level methylation matrices in one
command-line workflow. Its outputs integrate with common genomic tools for
quality control and downstream analysis.

## Highlights

- **Compact, cache-aware indexing.** Built with
  [libsais](https://github.com/IlyaGrebnov/libsais), bsbit encodes its combined
  three-letter search text in two bits per symbol and uses packed rank
  tables, sparse suffix-array samples, with bit-sliced counters. These
  cache-aware layouts reduce memory footprint and DRAM traffic, while bounded,
  in-place block merging avoids full-size intermediates and keeps index
  construction memory- and I/O-efficient.
- **Shared FM-index search.** One combined FM-index covers complementary
  sequence orientations, while dense seed lookup and
  [FMtree](https://doi.org/10.1093/bioinformatics/btx596)-inspired batched,
  multi-lane search share rank and locate work across reads.
- **Bit-parallel SIMD acceleration.** Myers bit-vector filters and
  SSE4.2/AVX2 SIMD kernels evaluate multiple candidate alignments per
  instruction, providing very high comparison throughput on modern x86-64
  processors.
- **Fast, accurate, and memory-efficient.** Candidate placements undergo
  complete edit-distance and CIGAR verification, while quality-aware evidence
  supports accurate placement and MAPQ. Under the same benchmark conditions,
  bsbit retained high accuracy and compact memory footprint while delivering
  about **1.5×** the throughput of the fastest competing aligner and
  **12–50×** that of the other tested aligners.
- **Beyond alignment.** bsbit provides integrated methylation quantification
  and SNV calling, including joint analysis and memory-efficient construction
  of cohort-scale methylation matrices.

## Installation and usage

For supported platforms, installation instructions, the quick start, workflow
guides, CLI options, input and output formats, and troubleshooting, see the
complete user guide:

[**wbvguo.github.io/bsbit/**](https://wbvguo.github.io/bsbit/)

## Citation

```text
Coming soon...
```

## License

bsbit is dual-licensed under the [MIT](LICENSE-MIT) or [Apache 2.0](LICENSE-APACHE) license.
Fixed third-party components retain their own licenses; see the
[third-party notices](external/licenses/THIRD_PARTY_NOTICES.md).
