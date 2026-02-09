---
template: home.html
title: Amplirust
---

<div class="feature-grid" markdown>

<div class="feature" markdown>
### :material-lightning-bolt: SIMD-Accelerated
Fast approximate primer matching powered by [sassy](https://github.com/RagnarGrootKoerkamp/sassy) with AVX2/SSE4/NEON support.
</div>

<div class="feature" markdown>
### :material-dna: IUPAC Support
Full IUPAC ambiguity code support (R, Y, S, W, K, M, B, D, H, V, N) in primer sequences.
</div>

<div class="feature" markdown>
### :material-file-multiple: Multi-Format
Reads FASTA and GenBank files with automatic format detection, including gzip/BGZF compressed inputs.
</div>

<div class="feature" markdown>
### :material-circle-outline: Circular Genomes
Handles plasmids and bacterial chromosomes with products that wrap around the origin.
</div>

<div class="feature" markdown>
### :material-cpu-64-bit: Multi-Threaded
Parallel file reading, searching, and compression for maximum throughput.
</div>

<div class="feature" markdown>
### :material-zip-box: Compression
Transparent gzip/BGZF support for both input and output with parallel decompression.
</div>

</div>

## Quick Start

```bash
# Install via conda
conda install bioconda::amplirust

# Run in-silico PCR
amplirust \
  -i genome.fasta \
  -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" \
  -o products.fasta
```

[:octicons-arrow-right-24: Getting Started](getting-started/installation.md){ .md-button }

## API Reference

Auto-generated API documentation from source code is available at the [API Reference](api/amplirust/).
