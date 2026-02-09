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

## IUPAC Ambiguity Codes

Primer sequences support standard IUPAC ambiguity codes:

| Code | Bases | Meaning |
|------|-------|---------|
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| S | G, C | Strong |
| W | A, T | Weak |
| K | G, T | Keto |
| M | A, C | Amino |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| H | A, C, T | Not G |
| V | A, C, G | Not T |
| N | A, C, G, T | Any |

These codes are expanded during matching, allowing a single primer to match multiple sequence variants.
