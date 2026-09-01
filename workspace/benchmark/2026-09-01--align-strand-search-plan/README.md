# Shared strand-search plan qualification — 2026-09-01

This evidence qualifies the responsibility refactor at `44f1d4c` on branch
`codex/align-strand-search-plan` against its parent `da5ea39`, plus the
source-surface cleanup at `920aa9a` against the qualified refactor.

The five-million-record GRCh38 fixture, frozen release binaries, and BAMs
remain outside Git under `/tmp/bsbit-align-strand-plan-20260901.FRNsk3`.
The checked-in files retain the commands, raw timings, internal work counts,
output digests, and interpretation needed to audit both stages.

Files:

- `report.md`: architecture decision, conclusions, and qualification scope.
- `results.tsv`: whole-program A/B timing and memory observations.
- `single-metrics.tsv`: invariant single-end classifications and work counts.
- `paired-metrics.tsv`: invariant paired-end classifications and work counts.
- `commands.md`: build, run, output-validation, and quality-gate commands.

The result is intentionally a maintainability qualification, not a speedup
claim: all four output modes are byte-identical, and matched user CPU changes
remain within about one percent.
