# Build index

Use `bsbit index` to build the complete alignment index from a reference genome
FASTA:

```bash
bsbit index \
  -r GRCh38.fa \
  -o GRCh38.bsbit \
  -t 8
```

## Parameters

<div class="cli-options" markdown>

| Option | Value | Default | Description |
|---|---|---|---|
| `-r`,<br>`--reference` | `PATH` | Required | Plain or BGZF-compressed reference genome FASTA used to build the index. |
| `-o`,<br>`--output` | `PATH` | Required | Path for the generated bsbit alignment index. |
| `-t`,<br>`--threads` | `N` | `1` | Number of threads used to build the index, from 1 to 64. |
| `-h`,<br>`--help` | — | None | Print help for `bsbit index` and exit. |

</div>

??? note "Parameter validation"

    The reference, output, and threads parameters each accept exactly one value
    in either short or long form. Unknown options, repeated options, a missing
    value, or a thread count outside the accepted range will be reported as
    errors.

??? tip "Faster alignment with a larger index"

    The default `balanced` index is recommended for most workflows. If the same
    reference index will be reused for many samples and alignment speed matters
    more than storage and memory use, select `fast` when building it:

    ```bash
    bsbit index \
      -r GRCh38.fa \
      -o GRCh38.fast.bsbit \
      -t 8 \
      --index-speed fast
    ```

    A `fast` index is larger and uses more memory during alignment, but lets
    bsbit locate candidate alignments more quickly. Alignment results are
    unchanged. In GRCh38 benchmarks with five million reads or read pairs,
    `fast` reduced alignment time by 6.5% for single-end and 10.4% for
    paired-end data, at a cost of about 1.44 GiB in index size and peak memory.

## Reference requirements

Plain [FASTA](../reference/file-formats.md#fasta-reference) is recommended. It
works for indexing without sidecars. For calling, an adjacent `.fai` is
recommended; without one, bsbit scans the FASTA once to build an in-memory
position table. BGZF-compressed FASTA is also accepted, but calling requires
both adjacent `.fai` and `.gzi` indexes. Create the indexes with `samtools
faidx`:

```bash
samtools faidx GRCh38.fa       # plain FASTA: creates .fai
samtools faidx GRCh38.fa.gz    # BGZF FASTA: creates .fai and .gzi
```

Ordinary gzip FASTA is not supported. Decompress it to plain FASTA or convert
it to BGZF:

```bash
gzip -cd GRCh38.fa.gz | bgzip -c > GRCh38.bgzf.fa.gz
samtools faidx GRCh38.bgzf.fa.gz
```

## Index output

`bsbit index` writes an opaque index bundle at `-o`. Pass the same path to
`bsbit align -x`. This index is used only for alignment; `bsbit call` reads the
original FASTA directly.

Treat the bundle as one artifact; downstream commands do not modify it after
generation. If an index already exists at the output path, bsbit replaces it
only after the new build succeeds. A failed build leaves the existing index
unchanged.

## Next

- [Align reads](alignment.md)
- [Input data](../reference/input-data.md)
- [CLI reference](../reference/cli.md#bsbit-index)
