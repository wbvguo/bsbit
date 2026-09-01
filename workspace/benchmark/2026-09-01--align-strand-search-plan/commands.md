# Reproduction commands

Commands were run from the `codex/align-strand-search-plan` worktree. Large
artifacts use `/tmp/bsbit-align-strand-plan-20260901.FRNsk3`.

## Freeze matched release binaries

Both the parent and candidate were built with:

```bash
CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit
```

The parent binary was copied before rebuilding the same target for the
candidate. Binary identity was recorded with:

```bash
sha256sum \
  /tmp/bsbit-align-strand-plan-20260901.FRNsk3/baseline-bsbit \
  /tmp/bsbit-align-strand-plan-20260901.FRNsk3/candidate-bsbit
```

## Five-million-record A/B runs

The common paths and CPU sets were:

```bash
EVIDENCE=/tmp/bsbit-align-strand-plan-20260901.FRNsk3
INDEX=/tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit
R1=/tmp/bsbit-current-benchmark-20260831/inputs/R1.fastq.gz
R2=/tmp/bsbit-current-benchmark-20260831/inputs/R2.fastq.gz
SINGLE_CPUS=0,2,4,6,8,10,12,14,16,18
PAIRED_CPUS=0,2,4,6,8,10,12,14,16,18,20,22,24,26
```

Each binary was run once per single-end mode and per non-directional paired
mode. Directional paired was repeated in reverse binary order:

```bash
taskset -c "$SINGLE_CPUS" "$EVIDENCE/candidate-bsbit" align \
  --index "$INDEX" --read1 "$R1" \
  --output-bam "$EVIDENCE/candidate-directional-single.bam" \
  --threads 8 --bam-threads 2 --bam-compression-level 1 --metrics

taskset -c "$SINGLE_CPUS" "$EVIDENCE/candidate-bsbit" align \
  --index "$INDEX" --read1 "$R1" --non-directional \
  --output-bam "$EVIDENCE/candidate-nondirectional-single.bam" \
  --threads 8 --bam-threads 2 --bam-compression-level 1 --metrics

taskset -c "$PAIRED_CPUS" "$EVIDENCE/candidate-bsbit" align \
  --index "$INDEX" --read1 "$R1" --read2 "$R2" \
  --output-bam "$EVIDENCE/candidate-directional-paired.bam" \
  --threads 10 --bam-threads 4 --bam-compression-level 1 --metrics

taskset -c "$PAIRED_CPUS" "$EVIDENCE/candidate-bsbit" align \
  --index "$INDEX" --read1 "$R1" --read2 "$R2" --non-directional \
  --output-bam "$EVIDENCE/candidate-nondirectional-paired.bam" \
  --threads 10 --bam-threads 4 --bam-compression-level 1 --metrics
```

The baseline commands replace `candidate` with `baseline`. `/usr/bin/time -v`
wrapped each invocation; the selected wall, user, system, and maximum-RSS
observations are preserved in `results.tsv`.

## Output validation

Run inside the raw evidence directory:

```bash
cmp -s baseline-directional-single.bam candidate-directional-single.bam
cmp -s baseline-nondirectional-single.bam candidate-nondirectional-single.bam
cmp -s baseline-directional-paired.bam candidate-directional-paired.bam
cmp -s baseline-nondirectional-paired.bam candidate-nondirectional-paired.bam
sha256sum *.bam
samtools quickcheck *.bam
```

All comparisons and quick checks succeeded.

## Quality gates

```bash
cargo check --locked --workspace --all-targets
cargo test --locked -p bsbit-align --all-targets
cargo test --locked -p bsbit-cli --all-targets
cargo test --locked --workspace --all-targets
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
cargo fmt --all -- --check
/mnt/conda/envs/dev/bin/python tests/tools/test_crate_boundaries.py
mkdocs build --strict
git diff --check
```
