# Installation

## Supported platforms

bsbit supports 64-bit Linux on x86-64 (Intel or AMD), either natively or under
WSL2, and on AArch64 (ARM64) with NEON.

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

Confirm that the CLI is available and inspect the installed version and
commands:

```bash
bsbit --version
bsbit --help
```

??? tip highlight "Verify runtime compatibility"
    CPU features visible inside WSL2, VMs, containers, and cluster jobs can
    differ from those available on the host. Check the environment where bsbit
    will actually run:

    ```bash
    bsbit cpu
    ```

    The command reports the detected features and selected backend. See the
    [`bsbit cpu` reference](../reference/cli.md#bsbit-cpu) for details.

<div class="next-step" markdown>

**Next:** Use the [quick start](quickstart.md) for common commands, or review
the [sequencing data support](workflow.md#sequencing-data-support).

</div>
