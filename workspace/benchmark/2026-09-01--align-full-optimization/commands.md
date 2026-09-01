# Reproduction commands

Commands were run from the `codex/align-full-optimization` worktree. Large
outputs use `/tmp/bsbit-align-full-optimization-20260901`.

## Production build

```bash
CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit
```

## Build a fast index

```bash
bsbit index \
  --reference /path/to/grch38.fa.gz \
  --output /path/to/grch38-sa8.bsbit \
  --threads 14 \
  --index-speed fast
```

Omitting `--index-speed` or using `balanced` emits the backward-compatible
stride-16 layout. The reader validates the format minor version, encoded
stride, and sampled-SA length before publishing the index.

## Full-corpus A/B harness

The immutable-run harness from the exploratory evidence was reused:

```bash
export BSBIT_PERF_RUN_ROOT=/tmp/bsbit-align-full-optimization-20260901/runs
export BSBIT_PERF_INDEX=/tmp/bsbit-align-full-optimization-20260901/sa8-index/grch38-sa8.bsbit

workspace/benchmark/2026-09-01--align-performance-experiments/run-case.sh \
  label /path/to/bsbit single 8 2 0,2,4,6,8,10,12,14,16,18

workspace/benchmark/2026-09-01--align-performance-experiments/run-case.sh \
  label /path/to/bsbit paired 10 4 \
  0,2,4,6,8,10,12,14,16,18,20,22,24,26
```

The final automatic paired invocation replaces explicit worker counts with:

```bash
bsbit align ... --total-threads 14
```

It selects 11 mapping + 3 BGZF workers for a balanced/stride-16 index and 10 +
4 for a fast/stride-8 index. Explicit `--threads` and `--bam-threads` remain
available and conflict with `--total-threads`.

## Final profile capture

```bash
CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt \
  -C force-frame-pointers=yes -C debuginfo=1 -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit

perf record -e cpu-clock:u -F 999 -g --call-graph fp -o single.perf.data -- \
  taskset -c 0,2,4,6,8,10,12,14,16,18 \
  /path/to/bsbit align --index /path/to/grch38-sa8.bsbit \
  --read1 /path/to/R1.fastq.gz --output-bam single.bam \
  --threads 8 --bam-threads 2

perf record -e cpu-clock:u -F 999 -g --call-graph fp -o paired.perf.data -- \
  taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26 \
  /path/to/bsbit align --index /path/to/grch38-sa8.bsbit \
  --read1 /path/to/R1.fastq.gz --read2 /path/to/R2.fastq.gz \
  --output-bam paired.bam --threads 10 --bam-threads 4 --metrics

perf report --stdio --percent-limit 0.25 --sort overhead,symbol \
  -i single.perf.data > single.flat.txt
```

## Final quality gates

```bash
cargo fmt --all -- --check
cargo test --locked -p bsbit-core --all-targets
cargo test --locked -p bsbit-align --all-targets
cargo test --locked -p bsbit-hts --all-targets
cargo test --locked -p bsbit-cli --all-targets
cargo test --locked -p bsbit-index --no-default-features \
  --features combined-index --all-targets
cargo test --locked -p bsbit-index --features index-construction --all-targets
python3 tests/tools/test_crate_boundaries.py
git diff --check
```
