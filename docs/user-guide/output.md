# Output Formats

Amplirust produces two output formats: FASTA sequences and an optional TSV statistics file.

## FASTA Output

Products are written with descriptive headers:

```
>reference_id:primer_name:1	pos=0-123	strand=+	len=124
ACGTACGTACGT...
>reference_id:primer_name_rc:2	pos=200-324	strand=-	len=125
TGCATGCATGCA...
```

### Header Format

```
{reference_id}:{primer_name}[_rc][_wrap]:{case_number}\tpos={start}-{end}\tstrand={+|-}\tlen={length}
```

| Field | Description |
|-------|-------------|
| `reference_id` | First whitespace-delimited token of the FASTA header, or LOCUS name for GenBank records (with fallback to accession, definition, or `unknown_N`) |
| `_rc` suffix | Product found on reverse complement strand |
| `_wrap` suffix | Product wraps around a circular genome |
| `case_number` | Increments per reference header (resets for each reference) |
| `strand` | `+` for forward, `-` for reverse complement match orientation |

Output sequences retain their strand orientation; the `strand` field indicates match orientation.

## TSV Statistics

The optional TSV output (`--tsv`) contains detailed information for each product:

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
| strand | `+` (forward) or `-` (reverse complement) |
| is_circular_wrap | `true` if product wraps around |
| product_seq | The extracted sequence |

## Worked Example

Given a genome with a 16S rRNA region and primers with 1 mismatch each:

```bash
amplirust -i genome.fasta \
  -p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT" \
  -k 2 -o products.fasta --tsv stats.tsv
```

**FASTA output:**

```
>genome_1:16S:1	pos=1024-2548	strand=+	len=1525
AGAGTTTGATCATGGCTCAG...AAGTCGTAACAAGGTAGCCGTA
```

**TSV row (selected columns):**

| Column | Value |
|--------|-------|
| amplicon_id | genome_1:16S:1 |
| primer_name | 16S |
| product_len | 1525 |
| fwd_start | 1024 |
| fwd_end | 1044 |
| fwd_mismatches | 0 |
| fwd_identity | 1.000 |
| rev_start | 2526 |
| rev_end | 2548 |
| rev_mismatches | 0 |
| rev_identity | 1.000 |
| strand | + |
