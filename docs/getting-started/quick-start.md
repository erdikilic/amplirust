# Quick Start

This guide walks you through your first in-silico PCR run with Amplirust.

## Prerequisites

Install Amplirust using one of the methods from the [Installation](installation.md) page:

```bash
conda install bioconda::amplirust
```

## Step 1: Prepare Input

You need a sequence file (FASTA or GenBank) and a primer pair. For this example, we'll search for 16S rRNA using universal primers.

Create a primer CSV file `primers.csv`:

```csv
name,forward,reverse
16S,AGAGTTTGATCMTGGCTCAG,TACGGYTACCTTGTTACGACTT
```

Or pass the primer directly on the command line:

```bash
-p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT"
```

## Step 2: Run Your First PCR

```bash
amplirust \
  -i genome.fasta \
  -p primers.csv \
  -o products.fasta \
  --tsv stats.tsv
```

## Step 3: Understand the Output

Amplirust produces two outputs:

**FASTA file** (`products.fasta`) -- extracted PCR product sequences:

```
>genome_1:16S:1	pos=100-1600	strand=+	len=1501
AGAGTTTGATCATGGCTCAG...
```

**TSV file** (`stats.tsv`) -- detailed match statistics including primer positions, edit distances, identity scores, and CIGAR strings.

**Summary** printed to stderr:

```
=== Amplirust Summary ===
Input sequences:  1
Primer pairs:     1
Products found:   1
```

## Next Steps

- [CLI Reference](../user-guide/usage.md) -- all command-line options
- [Primer Formats](../user-guide/primer-formats.md) -- inline and CSV primer input
- [Output Formats](../user-guide/output.md) -- FASTA header format and TSV columns
- [Practical Examples](../user-guide/examples.md) -- real-world workflows
