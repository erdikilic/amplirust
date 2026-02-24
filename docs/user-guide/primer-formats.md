# Primer Formats

Amplirust accepts primers either inline on the command line or from a CSV file.

## Command Line

Single primer pair:

```bash
-p "primer_name:FORWARD_SEQ:REVERSE_SEQ"
```

Multiple primer pairs (semicolon-separated):

```bash
-p "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT;ITS:TCCGTAGGTGAACCTGCGG:TCCTCCGCTTATTGATATGC"
```

## CSV File

Create a CSV file with the header `name,forward,reverse`:

```csv
name,forward,reverse
16S_V1V3,AGAGTTTGATCMTGGCTCAG,ATTACCGCGGCTGCTGG
16S_V3V4,CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC
ITS1,TCCGTAGGTGAACCTGCGG,GCTGCGTTCTTCATCGATGC
```

Then reference the file:

```bash
-p primers.csv
```

!!! tip "Which format to use?"
    Use inline primers for quick one-off searches. Use a CSV file when you have multiple primer pairs or want to reuse primers across runs.

## Pool Mode Primers

With `--pool`, primers are provided as individual sequences rather than forward/reverse pairs. Products are found between any two primers that match in correct orientation and distance.

### Inline Format

Single primer:

```bash
--pool -p "primer_name:SEQUENCE"
```

Multiple primers (semicolon-separated):

```bash
--pool -p "p1:AGAGTTTGATCMTGGCTCAG;p2:GWATTACCGCGGCKGCTG;p3:CCTACGGGNGGCWGCAG"
```

### CSV File (2 columns)

Create a CSV file with the header `name,sequence`:

```csv
name,sequence
27F,AGAGTTTGATCMTGGCTCAG
519R,GWATTACCGCGGCKGCTG
1492R,TACGGYTACCTTGTTACGACTT
```

Then reference the file:

```bash
--pool -p pool_primers.csv
```

!!! tip "When to use pool mode"
    Use pool mode (`--pool`) for multiplex panel design, MLST typing schemes, or exploratory primer screening where any two primers can form a product. Use `--pool-self-match` to also find products where the same primer acts as both forward and reverse.

## IUPAC Ambiguity Codes

Primer sequences support standard IUPAC ambiguity codes. These codes are expanded during matching, allowing a single primer to match multiple sequence variants.

| Code | Bases | Meaning |
|------|-------|---------|
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| S | G, C | Strong (3 H-bonds) |
| W | A, T | Weak (2 H-bonds) |
| K | G, T | Keto |
| M | A, C | Amino |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| H | A, C, T | Not G |
| V | A, C, G | Not T |
| N | A, C, G, T | Any |

### IUPAC Complement Pairs

When computing reverse complements, IUPAC codes follow these complement rules:

| Code | Complement | Reason |
|------|-----------|--------|
| R (A,G) | Y (C,T) | Purine ↔ Pyrimidine |
| S (G,C) | S (G,C) | Self-complementary |
| W (A,T) | W (A,T) | Self-complementary |
| K (G,T) | M (A,C) | Keto ↔ Amino |
| B (C,G,T) | V (A,C,G) | Not A ↔ Not T |
| D (A,G,T) | H (A,C,T) | Not C ↔ Not G |
| N | N | Any ↔ Any |
