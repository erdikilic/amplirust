# Introduction

Amplirust is a high-performance in-silico PCR tool for primer matching and product extraction from FASTA and GenBank sequences.

## Features

- **Fast approximate primer matching** using SIMD-accelerated algorithms (via [sassy](https://github.com/RagnarGrootKoerkamp/sassy))
- **IUPAC ambiguity code support** (R, Y, S, W, K, M, B, D, H, V, N) in primers
- **FASTA and GenBank input** with automatic format detection
- **Circular genome support** for plasmids and bacterial chromosomes
- **Reverse complement strand search** for comprehensive primer detection
- **Multi-threaded processing** for parallel file reading, searching, and compression
- **Gzip/BGZF support** for both input and output files (parallel decompression for BGZF)
- **Flexible primer input** via command line or CSV file

## Quick Start

```bash
# Install via conda
conda install bioconda::amplirust

# Run a simple PCR
amplirust -i genome.fasta -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" -o products.fasta
```

See the [Installation](./installation.md) chapter for all installation methods, or jump to [Usage](./usage.md) for detailed examples.

## API Reference

Auto-generated API documentation from source code is available at [API Reference](./api/amplirust/).
