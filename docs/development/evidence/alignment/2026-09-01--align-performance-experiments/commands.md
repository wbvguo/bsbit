# Reproduction commands

All commands were run from this worktree. Large outputs use
`/tmp/bsbit-align-perf-experiments-20260901`.

## Production build

The repository build wrapper rejected the baseline because its pinned root
`LICENSE-APACHE` checksum does not match clean commit `183d73e`; it was not
modified for this performance task. The equivalent production build was run
directly:

```bash
CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit
```

Each source state used a separate target directory or copied binary. Notable
preserved binaries are:

```text
/tmp/bsbit-align-perf-experiments-20260901/baseline-target/release/bsbit
/tmp/bsbit-align-perf-experiments-20260901/binaries/bounded-plus-paired-evaluator
/tmp/bsbit-align-perf-experiments-20260901/binaries/bounded-plus-single-fastpath
/tmp/bsbit-align-perf-experiments-20260901/binaries/fastpath-plus-packed-d3
/tmp/bsbit-align-perf-experiments-20260901/binaries/fastpath-sa8
```

## Full-corpus run

The checked-in harness refuses to overwrite a prior label. Example standard and
thread-mix invocations:

```bash
docs/development/evidence/alignment/2026-09-01--align-performance-experiments/run-case.sh \
  example-single /path/to/bsbit single 8 2 0,2,4,6,8,10,12,14,16,18

docs/development/evidence/alignment/2026-09-01--align-performance-experiments/run-case.sh \
  example-paired-11m3b /path/to/bsbit paired 11 3 \
  0,2,4,6,8,10,12,14,16,18,20,22,24,26
```

The stride-8 runs additionally set:

```bash
BSBIT_PERF_INDEX=/tmp/bsbit-align-perf-experiments-20260901/sa8-index/grch38-sa8.bsbit
```

## Correctness and repository checks

```bash
cargo fmt --all -- --check
cargo test --locked -p bsbit-align --all-targets
cargo test --locked -p bsbit-cli --all-targets
cargo test --locked -p bsbit-index --no-default-features \
  --features combined-index --all-targets
python3 tests/tools/test_crate_boundaries.py
git diff --check
```

Every full run also executes:

```bash
sha256sum output.bam
samtools quickcheck -v output.bam
samtools view -c output.bam
samtools view output.bam | awk '... mapped/unmapped/150M classification ...'
```

## Final profile build and capture

```bash
CARGO_INCREMENTAL=0 \
RUSTFLAGS='-C target-cpu=x86-64-v3 -C target-feature=+popcnt \
  -C force-frame-pointers=yes -C debuginfo=1 -A dead-code' \
CARGO_PROFILE_RELEASE_LTO=fat \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
cargo build --locked --release -p bsbit-cli --bin bsbit

perf record -e cpu-clock:u -F 999 -g --call-graph fp -o perf.data -- \
  taskset -c 0,2,4,6,8,10,12,14,16,18 \
  /path/to/bsbit align --index /path/to/stride16.bsbit \
  --read1 /path/to/R1.fastq.gz --output-bam output.bam \
  --threads 8 --bam-threads 2

perf report --stdio --percent-limit 0.5 --sort overhead,symbol \
  -i perf.data > flat.txt
perf report --stdio --children --percent-limit 0.5 --sort overhead,symbol \
  -i perf.data > callgraph.txt
```
