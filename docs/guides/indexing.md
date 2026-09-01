# Build index

Use one plain FASTA as the authoritative reference for alignment and calling.
Build its random-access index first, then build the complete `.bsbit`
alignment index:

```bash
samtools faidx GRCh38.fa

bsbit index \
  --reference GRCh38.fa \
  --output GRCh38.bsbit \
  --threads 8
```

Plain FASTA is the simplest format to share across the workflow. `bsbit index`
also accepts gzip or BGZF FASTA, but callers cannot randomly access ordinary
gzip FASTA; BGZF additionally needs a `.gzi`. Contig names must be unique, and
sequence is normalized to uppercase `A`, `C`, `G`, `T`, or `N`.

`--output` must be a new path and becomes the only opaque index handle supplied
to alignment. Physical construction and layout are internal. Alignment only
opens and validates the completed index; a
missing, stale, corrupt, or mismatched component is an error and is never
rebuilt during alignment. Downstream calling independently uses the indexed
original FASTA as its authoritative sequence.

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
- [Review index parameters](../reference/cli.md#bsbit-index)
