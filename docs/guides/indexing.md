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

| Option | Value | Default | Description |
|---|---|---|---|
| `-r`,<br>`--reference` | `PATH` | Required | Plain or BGZF-compressed reference genome FASTA used to build the index. |
| `-o`,<br>`--output` | `PATH` | Required | Path for the generated bsbit alignment index. |
| `-t`,<br>`--threads` | `N` | `1` | Number of threads used to build the index, from 1 to 64. |
| `-h`,<br>`--help` | — | None | Print help for `bsbit index` and exit. |

??? note "Parameter validation"

    The reference, output, and threads parameters each accept exactly one value
    in either short or long form. Unknown options, repeated options, a missing
    value, or a thread count outside the accepted range will be reported as
    errors.

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

`bsbit index` writes an opaque index bundle at `--output`. Pass the same path to
`bsbit align -i`. This index is used only for alignment; `bsbit call` reads the
original FASTA directly.

Treat the bundle as one artifact; downstream commands do not modify it after
generation. If an index already exists at the output path, bsbit replaces it
only after the new build succeeds. A failed build leaves the existing index
unchanged.

For repeated, throughput-oriented alignment where additional index storage and
mapping memory are acceptable, build the optional stride-8 index:

```bash
bsbit index \
  --reference GRCh38.fa \
  --output GRCh38.fast.bsbit \
  --threads 8 \
  --index-speed fast
```

The default `--index-speed balanced` keeps the stride-16 layout. Both modes
produce identical alignment decisions; `fast` stores roughly twice as many
sparse suffix-array samples so that locating a hit needs fewer LF steps.

## Next

- [Align reads](alignment.md)
- [Prepare input data](../reference/input-data.md)
- [Use the complete CLI reference](../reference/cli.md#bsbit-index)
