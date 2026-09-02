# `bsbit align` merged integration qualification — 2026-09-01

This directory is the durable qualification record for merging
`codex/non-directional-single-end` and `codex/single-sensitive-accuracy` into
`codex/align-optimized-integration`.

The code under test was merge commit `def3857`. The integration worktree was
already in the middle of merging the non-directional branch when consolidation
started, so that in-progress merge was completed safely as `b754ce7`; the
sensitive branch was then merged as `def3857`. Both source tips are ancestors
of the resulting integration history:

- `a237f51` -> `b754ce7` (non-directional single-end);
- `63d0f62` -> `def3857` (single-end sensitive accuracy).

Files:

- `report.md`: merge, correctness, BAM, and runtime conclusions;
- `results.tsv`: whole-program timing and output summary;
- `single-metrics.tsv`: exact built-in single-end metrics rows;
- `paired-metrics.tsv`: exact built-in paired-end metrics row;
- `commands.md`: exact build, run, and verification commands.

The two data rows in `single-metrics.tsv` are directional first and
non-directional second, matching the order in `results.tsv` and `report.md`.

The 5M fixtures, stride-8 index, BAM files, and release binary remain outside
Git under `/tmp`; their SHA-256 digests are recorded in `report.md`.
