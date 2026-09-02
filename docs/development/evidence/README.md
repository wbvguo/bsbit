# Compact evidence

This directory contains small, immutable evidence snapshots that are required
to support tracked documentation and therefore must be available in a clean
checkout.

`2026-08-29/` contains the adapter-recovery and sensitive-promotion tables
linked from the current performance evidence.

Complete run archives, logs, profiles, BAMs, and machine-local artifacts remain
under `workspace/runs/` or their originating
`workspace/worktree/YYYY-MM-DD--title/` attempt. Local cross-attempt synthesis
belongs under `workspace/notes/`. Product code, tests, and CI must not depend on
either ignored workspace path.
