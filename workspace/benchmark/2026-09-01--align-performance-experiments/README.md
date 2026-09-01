# `bsbit align` performance experiments — 2026-09-01

This directory records controlled experiments from branch
`codex/align-performance-experiments`, starting at commit `183d73e`.

Large BAM files, binaries, build trees, indexes, and `perf.data` files live under
`/tmp/bsbit-align-perf-experiments-20260901`. Small commands, timings, metrics,
checksums, and conclusions are retained here.

The primary full fixture is GRCh38 with five million 150 bp fragments:

- R1: `/tmp/bsbit-current-benchmark-20260831/inputs/R1.fastq.gz`
- R2: `/tmp/bsbit-current-benchmark-20260831/inputs/R2.fastq.gz`
- stride-16 index: `/tmp/bsbit-flattened-20260831/indices/bsbit/current.bsbit`

The similarly named index in the `bsbit-current-benchmark` tree is an older,
incompatible format and is deliberately not used.

Production binaries use `x86-64-v3,+popcnt`, fat LTO, and one codegen unit.
The standard comparison uses eight mapping workers, two BGZF workers, and the
ten physical CPUs `0,2,...,18`. See `report.md` for results and interpretation.

Files retained in git:

- `report.md`: conclusions, correctness checks, and experiment interpretation.
- `results.tsv`: selected full-corpus timing runs used by the report.
- `paired-metrics.tsv`: paired-end pipeline metrics for the thread/index study.
- `index-results.tsv`: index build and resident-memory costs.
- `commands.md`: reproducible build, run, validation, and profiling commands.
- `run-case.sh`: immutable-run harness that captures command, environment,
  `/usr/bin/time -v`, aligner metrics, BAM checksum, count, and classification.

Every run has its full evidence directory under
`/tmp/bsbit-align-perf-experiments-20260901/runs/<label>`. Runs marked invalid in
the tables were retained for audit but excluded from medians because the WSL
host showed an abnormal system-time or scheduling excursion.
