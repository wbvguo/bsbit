# Third-party distribution policy

This file describes the complete production scope matrix. Artifact assembly
generates a shorter `THIRD_PARTY_NOTICES.md` containing only the components in
the selected scope. Exact identities, source files, output paths, hashes, and
scope membership are recorded in `license-manifest.json`.

## Scope matrix

| Component group | `binary` | `source` | License material |
|---|---:|---:|---|
| HTSlib 1.24 | Yes | Yes | Complete HTSlib MIT/BSD text |
| htscodecs 1.6.7 | Yes | Yes | Complete BSD/public-domain/CC0 text |
| libsais 2.10.4 | Yes | Yes | Apache-2.0 text from the pinned source |
| Nine locked Rust registry packages | Yes | No | One shared root `LICENSE-APACHE` |
| Rust 1.94.0 standard library | Yes | No | Toolchain-generated `COPYRIGHT-library.html` |

HTSlib and htscodecs are both present in the `bsbit` binary: its linked static
archive contributes CRAM and codec objects even though the safe bsbit
API rejects CRAM input. Consequently both upstream license texts remain in the
binary scope.

The same `bsbit` binary owns index construction through `bsbit index`, so its
feature closure links the pinned libsais builder. The Apache-2.0 text therefore
belongs in both the binary and recursive-source scopes.

The locked Rust group is `block-buffer`, `cfg-if`, `cpufeatures`,
`crypto-common`, `digest`, `hybrid-array`, `libc`, `sha2`, and `typenum`.
They are redistributed under their Apache-2.0 option, share one unmodified
license text, and retain separate manifest rows only for versioned auditability.
Their audited registry sources contain no top-level `NOTICE` file.

Rust's standard library is statically present in the binary scope. Its
toolchain-generated copyright inventory covers the standard library and its
applicable exceptions, including code used for panic backtraces. Symbols such
as `miniz_oxide` and `adler2` in that runtime support do not implement bsbit's
FASTQ gzip path and are not packages in the workspace `Cargo.lock`.

Production FASTQ decoding uses HTSlib. Historical `zlib-rs`, `flate2`, and
`simd-adler32` experiments are not selectable by a production feature and are
absent from every release scope.

## Dynamically linked host components

The accepted x86_64 Linux `bsbit` binary uses host libraries such as system zlib,
libdeflate, liblzma, libbz2, libgomp, libm, libgcc_s, libc, and the ELF loader.
They are prerequisites rather than files in the standalone binary archive.
A package or container that redistributes them must inventory and satisfy the
licenses of the exact versions it includes.

## Project license

Project-owned work is available under `MIT OR Apache-2.0` with Copyright 2026
Wenbin Guo; see root `LICENSE-MIT` and `LICENSE-APACHE`. These project terms do
not relicense the components above.
