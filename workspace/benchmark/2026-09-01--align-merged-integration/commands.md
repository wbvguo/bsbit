# Reproduction commands

All commands were run from the `codex/align-optimized-integration` worktree.

## Merge topology

```bash
git show -s --format='%H %P %s' b754ce7
git show -s --format='%H %P %s' def3857
git merge-base --is-ancestor a237f51 def3857
git merge-base --is-ancestor 63d0f62 def3857
```

## Quality and release build

```bash
cargo check --workspace --all-targets
cargo test --locked --workspace --all-targets
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
cargo fmt --all -- --check
git diff --check
/mnt/conda/envs/dev/bin/python tests/tools/test_crate_boundaries.py
/mnt/conda/envs/dev/bin/python -m mkdocs build --strict \
  --site-dir /tmp/bsbit-docs-align-integration-merge

CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit
```

## Full-corpus runs

```bash
/usr/bin/time -v taskset -c 0,2,4,6,8,10,12,14,16,18 \
  target/release/bsbit align \
  --index /tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit \
  --read1 /tmp/bsbit-current-benchmark-20260831/inputs/R1.fastq.gz \
  --output-bam /tmp/bsbit-align-final-merge-20260901.ZgPvAd/directional-single.bam \
  --threads 8 --bam-threads 2 --bam-compression-level 1 --metrics

/usr/bin/time -v taskset -c 0,2,4,6,8,10,12,14,16,18 \
  target/release/bsbit align \
  --index /tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit \
  --read1 /tmp/bsbit-current-benchmark-20260831/inputs/R1.fastq.gz \
  --output-bam /tmp/bsbit-align-final-merge-20260901.ZgPvAd/non-directional-single.bam \
  --non-directional --threads 8 --bam-threads 2 \
  --bam-compression-level 1 --metrics

/usr/bin/time -v \
  taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26 \
  target/release/bsbit align \
  --index /tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit \
  --read1 /tmp/bsbit-current-benchmark-20260831/inputs/R1.fastq.gz \
  --read2 /tmp/bsbit-current-benchmark-20260831/inputs/R2.fastq.gz \
  --output-bam /tmp/bsbit-align-final-merge-20260901.ZgPvAd/paired.bam \
  --total-threads 14 --bam-compression-level 1 --metrics
```

## BAM verification

```bash
sha256sum /tmp/bsbit-align-final-merge-20260901.ZgPvAd/*.bam
samtools quickcheck -v /tmp/bsbit-align-final-merge-20260901.ZgPvAd/*.bam
samtools view -c BAM
samtools flagstat BAM
samtools view -H BAM
```
