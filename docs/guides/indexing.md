# Build index

Use `bsbit index` to build the complete alignment index from a reference genome
FASTA:

```bash
bsbit index \
  -r GRCh38.fa \
  -o GRCh38.bsbit \
  -t 8
```

## Common options

<div class="cli-options" markdown>

| Option | Argument | Default | Description |
| --- | --- | --- | --- |
| `-r`,<br>`--reference` | FASTA path | Required | Plain, gzip, or BGZF reference genome FASTA used to build the index |
| `-o`,<br>`--output` | Index path | Required | New path for the bsbit reference index |
| `-t`,<br>`--threads` | Positive integer | `1` | Number of indexing workers |

</div>

## Reference requirements

Uncompressed [FASTA](../reference/file-formats.md#fasta-reference) is
recommended for end-to-end workflows, although `bsbit index` can also read
gzip- and BGZF-compressed FASTA. A `.fai` is not required for the command but is
recommended to enable random access in `bsbit call`. Run the following command
for your reference format:

Uncompressed FASTA:

```bash
samtools faidx GRCh38.fa
```

BGZF-compressed FASTA (creates both `.fai` and `.gzi`):

```bash
samtools faidx GRCh38.fa.gz
```

??? note "Ordinary gzip versus BGZF"

    If `samtools faidx` reports `[E::fai_build3_core] Cannot index files
    compressed with gzip, please use bgzip`, the FASTA is compressed with
    standard gzip, which does not support indexed random access. Although
    `bsbit index` can read the file sequentially, we recommend converting it to
    BGZF and generating its `.fai` and `.gzi` indexes before running
    `bsbit index`, so the same reference is ready for `bsbit call`:

    ```bash
    gzip -cd downloaded-GRCh38.gz | bgzip -c > GRCh38.fa.gz
    samtools faidx GRCh38.fa.gz
    ```

## Index output

`bsbit index` writes the complete reference index at `-o`. Pass the same path
to `bsbit align -x`. This index is used only for alignment; `bsbit call` reads
the original FASTA directly.

Treat the bundle as one artifact; downstream commands do not modify it after
generation. If the reference FASTA changes, rerun `bsbit index` to generate a
new matching index before alignment.

## Next

- [Align reads](alignment.md)
- [Input data](../input-output/input-data.md)
- [CLI reference](../reference/cli.md#bsbit-index)
