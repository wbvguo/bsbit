# Installation

bsbit is built from source for 64-bit Linux on an Intel or AMD CPU with
x86-64-v3 support—roughly Intel Haswell or AMD Zen and newer. WSL2 provides
the supported Linux environment on x86-64 Windows; native Windows, macOS, and
ARM64 Linux are not qualified. Virtual machines must expose the required CPU
features.

??? tip "Check CPU compatibility"
    Run this inside the Linux, WSL2, VM, or container that will run bsbit:

    ```bash
    bsbit_cpu_result=$(
      /lib64/ld-linux-x86-64.so.2 --help 2>/dev/null \
        | grep -F 'x86-64-v3' || true
    )

    case "$bsbit_cpu_result" in
      *"(supported, searched)"*) echo "supported: x86-64-v3" ;;
      *x86-64-v3*) echo "unsupported: x86-64-v3 is not exposed" ;;
      *) echo "unknown: this glibc loader cannot report x86-64-v3" ;;
    esac
    ```

    `supported` means the current environment exposes the complete
    x86-64-v3 feature level. `unknown` usually means an older glibc or a
    different dynamic-loader path; it does not by itself mean that the CPU is
    unsupported.

## Install prerequisites

On Ubuntu or WSL2, install the native build and workflow tools:

```bash
sudo apt-get update
sudo apt-get install --yes \
  build-essential autoconf automake libtool pkg-config curl git \
  zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
  samtools tabix
```

Install Rust with [rustup](https://rustup.rs/) if `cargo` is unavailable:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. "$HOME/.cargo/env"
rustup default stable
```

Rust 1.89 or newer is required.

## Build

Clone the pinned native dependencies with the repository and build the complete
CLI:

```bash
git clone --recurse-submodules https://github.com/wbvguo/bsbit.git
cd bsbit
cargo build --locked --release -p bsbit-cli --bin bsbit

export PATH="$PWD/target/release:$PATH"
git rev-parse HEAD
```

The `PATH` change applies to the current shell. Repeat it in a new shell or
invoke `target/release/bsbit` directly.

For an existing clone with missing submodules, run:

```bash
git submodule update --init --recursive
```

## Verify

```bash
bsbit --version
bsbit --help
```

Keep the commit printed by `git rev-parse HEAD` with the analysis provenance.
bsbit is a source-first, pre-1.0 project and does not yet publish prebuilt
release artifacts or versioned documentation.

On WSL2, keep indexes and BAM files in the Linux filesystem rather than under
`/mnt/c/`. Use the [performance measurements](../performance-evidence.md) when
sizing storage and memory for a production reference.

Continue with the [quick start](quickstart.md), or review the
[supported workflows](workflow.md#supported-workflows) before using real data.
