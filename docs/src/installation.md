# Installation

## Conda (recommended)

The easiest way to install Amplirust is via [Bioconda](https://bioconda.github.io/):

```bash
conda install bioconda::amplirust
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

## FASTA Parser Selection

Amplirust supports two FASTA parsers via feature flags:

| Feature | Parser | Description |
|---------|--------|-------------|
| `parser_seqio` (default) | [seq_io](https://github.com/markschl/seq_io) | Fast, well-tested, minimal allocations |
| `parser_needletail` | [needletail](https://github.com/onecodex/needletail) | Very fast, used in production bioinformatics tools |

```bash
# Build with default parser (seq_io)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Build with needletail parser (potentially faster for large files)
RUSTFLAGS="-C target-cpu=native" cargo build --release --no-default-features --features parser_needletail

# Install with needletail
RUSTFLAGS="-C target-cpu=native" cargo install --path . --no-default-features --features parser_needletail
```
