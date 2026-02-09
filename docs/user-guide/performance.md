# Performance

## Build with Native CPU Optimizations

For best performance, compile with CPU-specific SIMD instructions:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This enables AVX2/SSE4 (x86) or NEON (ARM) acceleration for the approximate string matching engine.

## Thread Count

Adjust the number of threads based on your system. By default, Amplirust uses all available cores.

```bash
amplirust -t 8 -i large_genome.fasta -p primers.csv -o out.fasta
```

## Use BGZF Compression

For large single files, BGZF (blocked gzip) enables parallel decompression:

```bash
# Compress with bgzip (from htslib) instead of gzip
bgzip -c large_genome.fasta > large_genome.fasta.gz

# Amplirust auto-detects BGZF and uses parallel decompression
amplirust -i large_genome.fasta.gz -p primers.csv -o products.fasta.gz
```

Output `.gz` files are written in BGZF format (gzip-compatible).

## Increase Max Errors for Degenerate Primers

When working with degenerate primers or divergent sequences, increase the allowed edit distance:

```bash
amplirust -k 3 -i genome.fasta -p primers.csv -o out.fasta
```

!!! warning "Performance impact"
    Higher `-k` values increase search time. Use the minimum value that captures your expected matches, and combine with `--min-identity` for additional filtering without extra cost.

## Parser Selection

If processing very large FASTA files, the needletail parser may offer better performance:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release \
  --no-default-features --features parser_needletail
```

See [Installation](../getting-started/installation.md#fasta-parser-selection) for details on parser options.
