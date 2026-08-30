# Rust coverage-guided fuzz workspace

This independent Cargo workspace contains test-only libFuzzer targets. It is
not a product crate, is not a root-workspace member, and is never linked into
the aligner or release CLI.

## Pinned tools

- `cargo-fuzz 0.13.2`;
- `nightly-2026-07-31` (`rustc 1.99.0-nightly 8ab9fdff5`);
- `libfuzzer-sys 0.4.13` through the independent `tests/fuzz/Cargo.lock`;
- a C++17 compiler (the accepted host uses GCC/G++ 13.3.0).

The default stable Rust toolchain remains authoritative for product builds and
tests. The nightly toolchain is selected explicitly only for fuzz commands.
The runner suppresses only `dead_code` diagnostics to avoid a documented Rust
1.94 renderer crash on the accepted host's invalid WSL PTY dimensions.

## Targets

| Target | Input interpretation | Oracle |
|---|---|---|
| `fasta` | byte 0 selects `BufReader` capacity; remaining bytes are FASTA | repeated outcome equality, terminal replay after error, canonical semantic round trip |
| `fastq` | byte 0 selects capacity; remaining bytes are strict four-line FASTQ | repeated outcome equality, terminal replay after error, canonical semantic round trip |
| `paired_fastq` | two capacity bytes, little-endian split offset, then two concatenated FASTQ sources | repeated pair outcome equality, terminal replay, synchronized canonical round trip |

Text inputs are bounded to eight records, 512-byte physical lines, 384 bases
or qualities per record, and 1,024 total bases. The campaign also sets
libFuzzer `-max_len=4096`, a five-second per-input timeout, and a 2-GiB
sanitizer RSS limit.

## Running

For an ephemeral developer smoke:

```text
BSBIT_FUZZ_SECONDS_PER_TARGET=3 scripts/run-rust-fuzz.sh
```

For a source-closed formal capture, first commit the harness and use an absent
output directory:

```text
BSBIT_FUZZ_SECONDS_PER_TARGET=30 \
  scripts/run-rust-fuzz.sh artifacts/rust-coverage-fuzz
```

The formal runner requires a clean tracked worktree/index, runs offline, pins
the libFuzzer seed, validates all four logs, rejects sanitizer markers or crash
artifacts, checks pre/post source hashes and pinned submodule trees, and
publishes the artifact directory only after success.

LeakSanitizer is explicitly disabled on the accepted WSL/ptrace host because
that host does not provide a reliable leak-clean signal. AddressSanitizer
remains active for out-of-bounds, use-after-free, and related memory errors.
This campaign does not instrument HTSlib C code; the native ASan/UBSan mutation
harness covers that boundary separately.
