# Formal test support

Cross-crate, small, stable, redistributable fixtures live in `fixtures/`.
Most tests remain next to their owning crate under `crates/*/tests`; the
independent Rust coverage-guided workspace lives under `fuzz/`. Native-shim
fuzz targets remain with their owning shim under `crates/bsbit-hts/`.
Small fixture reproduction and result-evaluation utilities live under
`tools/`. Cross-component qualification and soak runners live under
`qualification/`. Test-only fuzz runners and their summarizers live beside the
independent workspace under `fuzz/`.

Nothing under `tests/` is invoked by a production binary at runtime. Build and
release entry points remain under `scripts/`.

Large datasets such as ERR2359938 FASTQ, reference genomes, BAM files, and
index images belong under `../workspace/datasets/`, never under `tests/`.
They must not be required by the default test suite or CI.
