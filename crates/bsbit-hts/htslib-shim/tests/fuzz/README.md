# Native HTSlib coverage targets

`native_fuzz.c` is compiled three times with
`BSBIT_NATIVE_FUZZ_MODE=1`, `2`, or `3` to exercise decoded input, SAM record,
and SAM header boundaries respectively. The campaign runner compiles HTSlib
1.24, the project shim, and this target with Clang sanitizer coverage; calling
an uninstrumented archive is explicitly insufficient.

Run an ephemeral campaign from the repository root with:

```text
BSBIT_NATIVE_FUZZ_SECONDS_PER_TARGET=1 \
  scripts/run-native-fuzz.sh
```

Pass one absent output directory for a committed-source formal capture. The
runner requires exact Clang/LLVM 18.1.3, builds from clean Git archives of the
pinned HTSlib and htscodecs revisions, applies strict warnings independently to
all three modes, and records archive/object sanitizer-coverage symbols before
executing any target.

Every input is capped at 4,096 bytes. The reader caps decoded output at 1 MiB,
checks terminal replay, and repeats the complete status/compression/byte-count
outcome. Record and header targets use only a private scratch output, validate
terminal writer semantics, and repeat status/final-size outcomes. Oracle
violations abort and therefore become ordinary libFuzzer failure artifacts.

The runner owns a private `TMPDIR` because the host's inherited Windows
`TMP`/`TEMP` paths may be read-only inside the WSL sandbox. AddressSanitizer and
UndefinedBehaviorSanitizer remain fail-fast; LeakSanitizer is disabled under
the documented WSL/ptrace host policy, so this campaign makes no leak-clean
claim.

The target owns no production API and is never linked into the Rust product.
All paths come from the runner-created `BSBIT_NATIVE_FUZZ_SCRATCH` directory.
