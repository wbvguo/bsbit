# Quick start

After [building bsbit](installation.md), run the bundled smoke test from the
repository root:

```bash
bash docs/examples/run-quickstart.sh
```

The script builds a tiny reference index, aligns four bisulfite read pairs,
checks the BAM, and verifies that all eight records have MAPQ 60. Its final line
prints the temporary output directory:

```text
quick-start smoke test passed: /tmp/...
```

The directory contains the reference and its FASTA index, `reference.bsbit`,
`alignment.bam`, and `alignment.summary.tsv`. This is a functional check of the
local build, not a performance benchmark or scientific validation dataset.

Next, read the [workflow guide](workflow.md) to choose an analysis path. The
[indexing](../guides/indexing.md) and [alignment](../guides/alignment.md) guides
show the commands to use with your own data.
