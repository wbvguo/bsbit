# bsbit HTSlib shim

This project-owned C library is the only intended C surface between the
private `sys` module in `bsbit-hts` and authenticated HTSlib 1.24. It is not a user
API and does not contain mapping, scoring, strand, tag, CLI, or publication
policy.

## Boundary

The public header exposes only opaque handles, integer status values, byte
buffers, byte counts, `errno`, and bounded ASCII diagnostics. No HTSlib struct,
native buffer, borrowed pointer, or integer sentinel crosses the ABI.

- `bsbit_hts_reader` decodes plain, gzip, and BGZF bytes through
  `bgzf_open`/`bgzf_read` and reports the content-derived compression class.
- `bsbit_hts_writer` parses the accepted canonical SAM header and one complete
  LF-terminated record at a time, then writes BAM through HTSlib.
- a failed read returns a zero byte count and makes later reads terminal;
- a failed record makes the writer terminal, so finalize cannot report success;
- close/finalize is explicit and checked; destroy is idempotent for null and
  performs only best-effort cleanup for an unclosed handle;
- handles are thread-confined. The ABI makes no `Send` or `Sync` promise.

The raw writer may truncate the path passed to HTSlib. Its only valid caller is
the future safe adapter, which must pass an exclusively created private staging
path and publish it through the accepted create-only lifecycle. Likewise, the
safe adapter must reject embedded-NUL paths, stdio, devices, URLs, plugin
schemes, and unsupported formats before calling this raw ABI.

## Version pin

Compilation requires `HTS_VERSION == 102400`, the CMake smoke requires an
exact `htslib=1.24` pkg-config record, and the runtime health check requires
`hts_version()` to equal `1.24`. Dependency source authentication is enforced
by the pinned recursive Git submodule revisions and
`scripts/check-native-sources.sh`. Release assembly copies the exact upstream
license from the submodule according to
`external/licenses/license-manifest.json`.

## Verification

`scripts/check-htslib-shim.sh NORMAL_PREFIX SANITIZED_SOURCE` performs:

1. a strict C11 CMake build with project warnings promoted to errors;
2. an 84-assertion smoke covering ABI health, misleading suffixes, plain/gzip/
   BGZF equality, EOF, terminal truncation errors, close errors, canonical
   SAM-to-BAM encoding, independent BAM read-back, and writer failure state;
3. a deterministic mutation smoke over 512 truncated/bit-mutated gzip inputs,
   512 truncated/bit-mutated canonical SAM records, and 256 similarly mutated
   SAM headers; every case is bounded and checks the public terminal-state and
   cleanup contract;
4. a Linux link-time `--wrap` fault smoke with 51 assertions over 13 allocation,
   open, classify, read, close, parse, write, and finalize failpoints; wrappers
   call the real close before returning an injected close failure, so ownership
   counts are checked without manufacturing a leak;
5. an exported-symbol allow-list check; mutation/fault controls are test-only
   symbols and never enter the production library;
6. all three smokes with the shim and an instrumented static HTSlib under ASan and
   UBSan (`detect_leaks=0` is explicit because the reference host is under
   ptrace);
7. the pinned native-submodule source check.

Both build trees and all artifacts are created in a private temporary directory.
No generated file is written below this source directory.
