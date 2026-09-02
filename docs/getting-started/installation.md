# Installation

## Supported platforms

bsbit supports 64-bit Linux, either natively or inside WSL2, on an Intel or AMD
CPU with x86-64-v3 support—roughly Intel Haswell or AMD Zen and newer. Virtual
machines must expose the required CPU features.

??? tip highlight "Check CPU compatibility"
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

## Install from Bioconda

!!! info "Coming soon"
    Bioconda installation is coming soon. It will be the recommended
    installation method, especially on shared systems and HPC clusters.

## Install from source

### Requirements

The following tools and libraries are required to build bsbit and run the
documented workflow:

- Rust 1.89 or later
- A C/C++ toolchain with OpenMP support
- Autoconf, Automake, Libtool, Make, pkg-config, Git, and curl
- Development libraries for zlib, bzip2, liblzma, and libdeflate
- samtools and tabix

Choose one of the following ways to prepare a complete build environment.

??? tip highlight "Using Conda or Mamba"
    Conda can provide the compilers, Rust toolchain, native libraries, and
    workflow tools in an isolated user-space environment:

    ```bash
    conda create -n bsenv -c conda-forge -c bioconda \
      'rust>=1.89' c-compiler cxx-compiler make \
      autoconf automake libtool pkg-config git curl \
      zlib bzip2 xz libdeflate samtools htslib
    conda activate bsenv
    ```

    `mamba` can be used in place of `conda`. Activate `bsenv` again before
    building or running bsbit in a new shell.

??? tip highlight "Using system tools"
    On Ubuntu or WSL2, install the native dependencies and workflow tools with
    APT:

    ```bash
    sudo apt-get update
    sudo apt-get install --yes \
      build-essential autoconf automake libtool pkg-config curl git \
      zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
      samtools tabix
    ```

    If `cargo` is unavailable, install Rust with
    [rustup](https://rustup.rs/):

    ```bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    . "$HOME/.cargo/env"
    rustup default stable
    ```

### Build and install bsbit

Clone the pinned native dependencies and install the complete CLI using either
the Conda or system environment:

```bash
git clone --recurse-submodules https://github.com/wbvguo/bsbit.git
cd bsbit
cargo install --locked --path crates/bsbit-cli \
  --root "${CONDA_PREFIX:-$HOME/.cargo}"
```

After installation, the `bsbit` executable should be available on `PATH`.

### Verify the installation

```bash
bsbit --version
```

<div class="next-step" markdown>

**Next:** Use the [quick start](quickstart.md) for common commands, or review
the [sequencing data support](workflow.md#sequencing-data-support).

</div>
