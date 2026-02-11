# Amplirust

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://erdikilic.github.io/amplirust/)
[![codecov](https://codecov.io/gh/erdikilic/amplirust/branch/quality-overhaul/graph/badge.svg?token=7R102309W1)](https://codecov.io/gh/erdikilic/amplirust)
[![Conda Version](https://anaconda.org/bioconda/amplirust/badges/version.svg)](https://anaconda.org/bioconda/amplirust)
[![Conda Downloads](https://anaconda.org/bioconda/amplirust/badges/downloads.svg)](https://anaconda.org/bioconda/amplirust)
[![Conda Platforms](https://anaconda.org/bioconda/amplirust/badges/platforms.svg)](https://anaconda.org/bioconda/amplirust)

A high-performance in-silico PCR tool for primer matching and product extraction from FASTA and GenBank sequences.

## Features

- **Fast approximate primer matching** using SIMD-accelerated algorithms (via [sassy](https://github.com/RagnarGrootKoerkamp/sassy))
- **IUPAC ambiguity code support** (R, Y, S, W, K, M, B, D, H, V, N) in primers
- **FASTA and GenBank input** with automatic format detection
- **Circular genome support** for plasmids and bacterial chromosomes
- **Reverse complement strand search** for comprehensive primer detection
- **Multi-threaded processing** for parallel file reading, searching, and compression
- **Gzip/BGZF support** for both input and output files (parallel decompression for BGZF)
- **Flexible primer input** via command line or CSV file

## Installation

### Conda (recommended)

```bash
conda install bioconda::amplirust
```

### From Source

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

### FASTA Parser Selection

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

## Usage

```
amplirust [OPTIONS] --input <INPUT> --primers <PRIMERS>
```

### Basic Examples

```bash
# Simple PCR with inline primer
amplirust -i genome.fasta -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" -o products.fasta

# Multiple input files with glob pattern
amplirust -i "genomes/*.fna.gz" -p primers.csv -o results.fasta

# Circular genome (e.g., plasmid)
amplirust -i plasmid.fasta -p "ori:ACGTACGT:TGCATGCA" --circular -o products.fasta

# Search both strands with detailed output
amplirust -i genome.fasta -p primers.csv --search-rc -o products.fasta --tsv stats.tsv

# Gzip compressed output
amplirust -i genome.fasta.gz -p primers.csv -o products.fasta.gz

# GenBank input
amplirust -i sequence.gbk -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" -o products.fasta

# Mixed formats via glob (FASTA and GenBank files auto-detected)
amplirust -i "data/*" -p primers.csv -o products.fasta
```

### Options

#### Input Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input <FILES>` | | Input sequence files (comma-separated, glob patterns supported). Supports FASTA (`.fasta`, `.fa`, `.fna`, `.ffn`, `.fas`) and GenBank (`.gb`, `.gbk`, `.gbff`, `.genbank`, `.gbf`), including gzip-compressed (`.gz`). Unrecognized file types in glob results are skipped with a warning. |
| `-p, --primers <PRIMERS>` | | Primers as `name:forward:reverse` or path to CSV file |
| `--circular` | false | Treat sequences as circular genomes |
| `--max-decompression-size <N>` | 4294967296 | Maximum decompressed file size in bytes (0 = unlimited) |

#### Matching Options

| Option | Default | Description |
|--------|---------|-------------|
| `-k, --max-errors <N>` | 2 | Maximum edit distance for primer matching |
| `--min-identity <FLOAT>` | 0 (disabled) | Minimum identity threshold (0.0-1.0) |
| `--search-rc` | false | Also search reverse complement strand |
| `-t, --threads <N>` | auto | Number of threads (0 = auto-detect) |

#### Product Options

| Option | Default | Description |
|--------|---------|-------------|
| `--min-len <N>` | 50 | Minimum PCR product length |
| `--max-len <N>` | 5000 | Maximum PCR product length |
| `--trim-primers` | false | Remove primer sequences from output |
| `--max-n-fraction <F>` | 0.1 | Maximum fraction of N bases allowed in product sequence (0.0 - 1.0) |

#### Output Options

| Option | Description |
|--------|-------------|
| `-o, --output <FILE>` | Output FASTA file (.gz for compression) |
| `--tsv <FILE>` | TSV file with detailed statistics |
| `-v, --verbose` | Increase verbosity (-v, -vv, -vvv) |
| `-q, --quiet` | Suppress progress bar output |
| `--remove-duplicates` | Remove duplicate product sequences per reference (canonicalized by reverse complement) |

## Progress Bar and Summary

Amplirust shows a progress bar while reading multiple input files (suppressed with `--quiet` or any `-v` flag). It always prints a run summary to stderr after processing.

```
⠋ [00:00:02] [########################################] 15/15 files (0s)

=== Amplirust Summary ===
Input sequences:  1250
Primer pairs:     3
Products found:   47

Products by primer:
  16S: 45
  ITS: 2

Products by reference:
  genome_1: 30
  genome_2: 17

Output written to: products.fasta
TSV written to: stats.tsv
```

Use `--quiet` to suppress the progress bar (useful for scripting). Progress bars are also disabled when any verbosity flag is set (`-v`, `-vv`, `-vvv`).

## Primer Input Formats

### Command Line

Single primer pair:
```bash
-p "primer_name:FORWARD_SEQ:REVERSE_SEQ"
```

Multiple primer pairs (semicolon-separated):
```bash
-p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT;ITS:TCCGTAGGTGAACCTGCGG:TCCTCCGCTTATTGATATGC"
```

### CSV File

Create a CSV file with header:
```csv
name,forward,reverse
16S_V1V3,AGAGTTTGATCMTGGCTCAG,ATTACCGCGGCTGCTGG
16S_V3V4,CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC
ITS1,TCCGTAGGTGAACCTGCGG,GCTGCGTTCTTCATCGATGC
```

Then use:
```bash
-p primers.csv
```

## Output Formats

### FASTA Output

Products are written with descriptive headers:
```
>reference_id:primer_name:1	pos=0-123	strand=+	len=124
ACGTACGTACGT...
>reference_id:primer_name_rc:2	pos=200-324	strand=-	len=125
TGCATGCATGCA...
```

Header format: `{reference_id}:{primer_name}[_rc][_wrap]:{case_number}\tpos={start}-{end}\tstrand={+|-}\tlen={length}`
- `reference_id` is the first whitespace-delimited token of the FASTA header, or the LOCUS name for GenBank records (with fallback to accession, definition, or `unknown_N`)
- `_rc` suffix indicates product from reverse complement strand
- `_wrap` suffix indicates product wraps around circular genome
- `case_number` increments per reference header (resets for each reference)
- Output sequences retain their strand orientation; `strand` indicates match orientation.

### TSV Statistics

The TSV output contains detailed information for each product:

| Column | Description |
|--------|-------------|
| amplicon_id | Header using `reference_id` with case number |
| reference_id | Original sequence header |
| source_file | Input file path |
| primer_name | Primer pair name |
| product_len | Output sequence length |
| full_len | Full product length (before trimming) |
| fwd_start | Forward primer match start (0-based) |
| fwd_end | Forward primer match end |
| fwd_mismatches | Edit distance for forward primer |
| fwd_identity | Identity percentage for forward |
| fwd_cigar | CIGAR string for forward alignment |
| rev_start | Reverse primer match start |
| rev_end | Reverse primer match end |
| rev_mismatches | Edit distance for reverse primer |
| rev_identity | Identity percentage for reverse |
| rev_cigar | CIGAR string for reverse alignment |
| strand | + (forward) or - (reverse complement) |
| is_circular_wrap | true if product wraps around |
| product_seq | The extracted sequence |

## Performance Tips

1. **Use native CPU optimizations** for best performance:
   ```bash
   RUSTFLAGS="-C target-cpu=native" cargo build --release
   ```

2. **Adjust thread count** based on your system:
   ```bash
   amplirust -t 8 -i large_genome.fasta -p primers.csv -o out.fasta
   ```

3. **Use BGZF compression** for large single files (enables parallel decompression):
   ```bash
   # Compress with bgzip (from htslib) instead of gzip
   bgzip -c large_genome.fasta > large_genome.fasta.gz

   # Amplirust auto-detects BGZF and uses parallel decompression
   amplirust -i large_genome.fasta.gz -p primers.csv -o products.fasta.gz  # input: bgzf
   ```
   Output `.gz` files are written in BGZF format (gzip-compatible).

4. **Increase max-errors** for degenerate primers or divergent sequences:
   ```bash
   amplirust -k 3 -i genome.fasta -p primers.csv -o out.fasta
   ```

## Examples

### Extract 16S rRNA Regions

```bash
amplirust \
  -i bacterial_genomes/*.fasta \
  -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" \
  -k 2 \
  --min-len 1400 \
  --max-len 1600 \
  -o 16s_sequences.fasta \
  --tsv 16s_stats.tsv \
  -v
```

### Search Plasmid with Circular Mode

```bash
amplirust \
  -i plasmid.fasta \
  -p primers.csv \
  --circular \
  --search-rc \
  -o plasmid_products.fasta \
  -vv
```

### High-Sensitivity Search

```bash
amplirust \
  -i divergent_genome.fasta \
  -p primers.csv \
  -k 4 \
  --min-identity 0.75 \
  --search-rc \
  -o products.fasta
```

## License

MIT License
