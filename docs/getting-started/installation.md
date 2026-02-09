# Installation

## Conda (recommended)

The easiest way to install Amplirust is via [Bioconda](https://bioconda.github.io/):

=== "Conda"

    ```bash
    conda install bioconda::amplirust
    ```

=== "Mamba"

    ```bash
    mamba install bioconda::amplirust
    ```

Pre-built binaries are available for Linux (x86_64, aarch64) and macOS (x86_64, aarch64).

### Verify installation

```bash
amplirust --version
```

## From Source

Requires Rust 1.85+ (2024 edition).

```bash
# Clone the repository
git clone https://github.com/erdikilic/amplirust.git
cd amplirust

# Build with SIMD optimizations (recommended)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# The binary will be at target/release/amplirust
```

### Quick Install

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --path .
```

!!! tip "Native CPU optimizations"
    Always use `RUSTFLAGS="-C target-cpu=native"` when building from source. This enables AVX2/SSE4 (x86) or NEON (ARM) acceleration for the approximate string matching engine, providing significant performance improvements.

## FASTA Parser Selection

Amplirust supports two FASTA parsers via feature flags:

| Feature | Parser | Description |
|---------|--------|-------------|
| `parser_seqio` (default) | [seq_io](https://github.com/markschl/seq_io) | Fast, well-tested, minimal allocations |
| `parser_needletail` | [needletail](https://github.com/onecodex/needletail) | Very fast, used in production bioinformatics tools |

=== "Default (seq_io)"

    ```bash
    RUSTFLAGS="-C target-cpu=native" cargo build --release
    ```

=== "Needletail"

    ```bash
    RUSTFLAGS="-C target-cpu=native" cargo build --release \
      --no-default-features --features parser_needletail
    ```

=== "Install with needletail"

    ```bash
    RUSTFLAGS="-C target-cpu=native" cargo install --path . \
      --no-default-features --features parser_needletail
    ```
