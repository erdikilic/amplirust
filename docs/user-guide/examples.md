# Practical Examples

## Extract 16S rRNA Regions

Search bacterial genomes for 16S rRNA using universal primers with expected product length constraints:

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

## Search Plasmid with Circular Mode

For plasmids and circular chromosomes, enable `--circular` to find products that wrap around the origin:

```bash
amplirust \
  -i plasmid.fasta \
  -p primers.csv \
  --circular \
  --search-rc \
  -o plasmid_products.fasta \
  -vv
```

## High-Sensitivity Search

Increase error tolerance and lower identity threshold for divergent sequences:

```bash
amplirust \
  -i divergent_genome.fasta \
  -p primers.csv \
  -k 4 \
  --min-identity 0.75 \
  --search-rc \
  -o products.fasta
```

## Batch Processing with Globs

Process all FASTA files in a directory (including gzip-compressed):

```bash
amplirust \
  -i "genomes/*.fna.gz" \
  -p primers.csv \
  -o all_products.fasta \
  --tsv all_stats.tsv
```

Process multiple directories:

```bash
amplirust \
  -i "project_a/*.fasta,project_b/*.fna" \
  -p primers.csv \
  -o combined_products.fasta
```

## Compressed I/O Workflows

Amplirust transparently handles gzip and BGZF compressed files:

```bash
# Compressed input and output
amplirust \
  -i genomes.fasta.gz \
  -p primers.csv \
  -o products.fasta.gz

# BGZF input for parallel decompression (large files)
bgzip -c large_genome.fasta > large_genome.fasta.gz
amplirust -i large_genome.fasta.gz -p primers.csv -o products.fasta.gz
```

Output `.gz` files are written in BGZF format (gzip-compatible).

## GenBank Input Processing

Amplirust auto-detects GenBank format from file extensions (`.gb`, `.gbk`, `.gbff`, `.genbank`, `.gbf`):

```bash
amplirust \
  -i sequence.gbk \
  -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" \
  -o products.fasta
```

Mixed formats are also supported via glob patterns -- unrecognized file types are skipped with a warning:

```bash
amplirust -i "data/*" -p primers.csv -o products.fasta
```

## Primer Pool Mode (Multiplex Screening)

Find products between any combination of primers from a pool:

```bash
# Basic pool mode with inline primers
amplirust --pool \
  -i genome.fasta \
  -p "27F:AGAGTTTGATCMTGGCTCAG;519R:GWATTACCGCGGCKGCTG;1492R:TACGGYTACCTTGTTACGACTT" \
  -o products.fasta
```

Pool mode with a CSV file:

```bash
# CSV with 2 columns: name,sequence
amplirust --pool \
  -i bacterial_genomes/*.fasta \
  -p pool_primers.csv \
  -k 2 \
  -o pool_products.fasta \
  --tsv pool_stats.tsv
```

Allow self-matching (same primer as both forward and reverse):

```bash
amplirust --pool --pool-self-match \
  -i genome.fasta \
  -p "IR:CCTGCAGGCATGCAAGCTT" \
  -o self_products.fasta
```

Pool mode combined with circular genome search:

```bash
amplirust --pool \
  -i plasmid.fasta \
  -p pool_primers.csv \
  --circular \
  --search-rc \
  -o plasmid_pool_products.fasta
```

## TSV Analysis Workflow

Export detailed statistics and analyze with standard tools:

```bash
# Run PCR with TSV output
amplirust \
  -i genomes/*.fasta \
  -p primers.csv \
  -o products.fasta \
  --tsv stats.tsv

# Filter high-identity hits
awk -F'\t' 'NR==1 || ($10 >= 0.95 && $15 >= 0.95)' stats.tsv > high_identity.tsv

# Count products per primer
awk -F'\t' 'NR>1 {print $4}' stats.tsv | sort | uniq -c | sort -rn
```
