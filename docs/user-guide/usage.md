# CLI Reference

```
amplirust [OPTIONS] --input <INPUT> --primers <PRIMERS>
```

## Input Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input <FILES>` | | Input sequence files (comma-separated, glob patterns supported). Supports FASTA (`.fasta`, `.fa`, `.fna`, `.ffn`, `.fas`) and GenBank (`.gb`, `.gbk`, `.gbff`, `.genbank`, `.gbf`), including gzip-compressed (`.gz`). Unrecognized file types in glob results are skipped with a warning. |
| `-p, --primers <PRIMERS>` | | Primers as `name:forward:reverse` or path to CSV file |
| `--circular` | false | Treat sequences as circular genomes |
| `--max-decompression-size <N>` | 4294967296 | Maximum decompressed file size in bytes (0 = unlimited) |

## Matching Options

| Option | Default | Description |
|--------|---------|-------------|
| `-k, --max-errors <N>` | 2 | Maximum edit distance for primer matching |
| `--min-identity <FLOAT>` | 0 (disabled) | Minimum identity threshold (0.0-1.0) |
| `--search-rc` | false | Also search reverse complement strand |
| `-t, --threads <N>` | auto | Number of threads (0 = auto-detect) |

## Product Options

| Option | Default | Description |
|--------|---------|-------------|
| `--min-len <N>` | 50 | Minimum PCR product length |
| `--max-len <N>` | 5000 | Maximum PCR product length |
| `--trim-primers` | false | Remove primer sequences from output |
| `--max-n-fraction <F>` | 0.1 | Maximum fraction of N bases allowed in product sequence (0.0-1.0) |

## Output Options

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
