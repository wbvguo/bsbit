# bsbit

`bsbit` is an ultrafast, memory-efficient toolkit for bisulfite-sequencing
alignment, methylation and SNP calling, and cohort-matrix construction. The
qualified release target is 64-bit Linux on an Intel or AMD CPU with
x86-64-v3 support—roughly Intel Haswell or AMD Zen and newer. It builds a
memory-efficient reference index, reads plain or gzip FASTQ directly, and
publishes create-only BAM through HTSlib.

Start with the [documentation site](docs/index.md). It is the authoritative
source for installation, workflows, command arguments, output contracts,
troubleshooting, and current qualification evidence. Historical experiments do
not compile into the product; their lifecycle is documented under
[development](docs/development/feature-lifecycle.md).

## Prerequisites

Building requires Rust 1.89 or newer, a C toolchain, Autotools, Make,
pkg-config, OpenMP, and development libraries for zlib, bzip2, liblzma, and
libdeflate. On Debian or Ubuntu:

```sh
sudo apt-get install build-essential autoconf automake libtool pkg-config \
  zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev
```

HTSlib/htscodecs and libsais are pinned Git submodules under `external/`.
Clone with `--recurse-submodules`, or initialize an existing checkout with:

```sh
git submodule update --init --recursive
```

See [installation](docs/getting-started/installation.md) for the supported
platform contract and complete setup notes.

## Build

Build the complete indexing, alignment, calling, and combining command surface:

```sh
cargo build --locked --release -p bsbit-cli --bin bsbit
```

For an audited x86-64-v3 fat-LTO `bsbit` binary, use
`scripts/build-bsbit.sh`. It writes build products and provenance below
ignored `build/` unless given an explicit temporary output path.

## Minimal index and alignment flow

Create the complete reference index:

```sh
target/release/bsbit index \
  --reference GRCh38.fa.gz \
  --output GRCh38.bsbit \
  --threads 8
```

Use the same command for either read layout. This paired-end example supplies
both mates; alignment only opens and validates the opaque index created above
and never builds or changes it:

```sh
target/release/bsbit align \
  --index GRCh38.bsbit \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output-bam sample.bam \
  --threads 8 \
  --bam-threads 2
```

For directional single-end input, supply only `--read1` (or `-1`) and no
paired-only options. The single-end path writes numeric MAPQ from the retained
search evidence and declares caller-compatible provenance; see the
[single-end alignment guide](docs/guides/alignment.md).

Omit a mode flag for the default mode or add `--sensitive` for the qualified
maximum-recall mode. Search, rescue, reporting, and MAPQ policies are fixed by
each mode rather than assembled from experimental switches. Use `--help` for
the exact command surface, and see [indexing](docs/guides/indexing.md),
[alignment](docs/guides/alignment.md) before changing resource or output
settings.

## Outputs and downstream tools

The default `minimal` BAM contract emits standard BAM core fields plus literal
`NM` and bisulfite-strand `XG`. Select `--output-contract bismark` only when a
consumer requires Bismark-compatible `MD`, `XM`, `XR`, and `XG` tags. The
default BAM retains mapped MAPQ-0 representatives and unmapped primary records;
`--mapped-only` removes only truly unmapped records.

`bsbit call meth`, `bsbit call snp`, and `bsbit call joint` consume a
coordinate-sorted, duplicate-marked, indexed BAM plus its indexed reference.
`bsbit combine` merges named extended bedMethyl files into count and/or level
matrices. Output schemas and preprocessing contracts are maintained in:

- [workflow outputs](docs/outputs/index.md)
- [methylation calls](docs/outputs/methylation.md)
- [variant and joint calls](docs/outputs/variant-calling.md)
- [methylation matrices](docs/outputs/methylation-matrices.md)

The full command reference is [docs/reference/cli.md](docs/reference/cli.md).

## Validation

Run the normal source checks from the repository root:

```sh
scripts/check-native-sources.sh
cargo fmt --all -- --check
cargo test --locked --workspace
```

Formal fuzz, native-boundary, platform-publication, and release-soak entry
points are indexed in [scripts/README.md](scripts/README.md). Large references,
FASTQ inputs, BAM outputs, profiling data, and benchmark runs remain under
ignored `workspace/` or `agent/`; they are not required by the default suite.

## Project boundaries and evidence

Tracked product code lives in `crates/`, durable small fixtures in `tests/fixtures/`,
and current user/developer contracts in `docs/`. See
[repository layout](docs/repository-layout.md) and
[architecture](docs/architecture.md) for crate ownership and dependency
direction.

Performance numbers, support boundaries, MAPQ interpretation, and known
differences change with qualification evidence and therefore live only in the
[current performance evidence](docs/performance-evidence.md),
[supported workflows](docs/getting-started/workflow.md#supported-workflows), and
[known differences](docs/known-differences.md) pages.

## Build the documentation

Install the locked documentation dependency group in a virtual environment or
use the existing Conda environment, then preview the site on port 8800:

```sh
python3 -m venv .venv-docs
source .venv-docs/bin/activate
python -m pip install --upgrade "pip>=25.1"
python -m pip install --group docs
mkdocs serve
```

Open <http://127.0.0.1:8800/bsbit/>. Set `DEV_ADDR` only when another listen
address is required. Before publishing, run `mkdocs build --strict`; CI runs
the same strict build.

Project-owned code is dual-licensed under MIT or Apache-2.0. Fixed third-party
components retain their own licenses; see
[`external/licenses/THIRD_PARTY_NOTICES.md`](external/licenses/THIRD_PARTY_NOTICES.md).
