# Algorithm

This page describes how Amplirust performs in-silico PCR.

## Overview

In-silico PCR simulates the polymerase chain reaction computationally. Given a pair of primers (forward and reverse) and a target sequence, the algorithm:

1. Searches for approximate matches of the forward primer in the sequence
2. For each forward match, searches for the reverse complement of the reverse primer downstream
3. Extracts the region between (and including) the primer binding sites as a PCR product
4. Filters products by length, identity, and N-base content

## SIMD-Accelerated Approximate Matching

Amplirust uses the [sassy](https://github.com/RagnarGrootKoerkamp/sassy) library for approximate string matching. Sassy provides SIMD-accelerated (AVX2, SSE4, NEON) pattern matching with edit distance constraints, making it significantly faster than naive approaches.

The search finds all occurrences of a primer in the target sequence within a specified maximum edit distance (`-k`). Edit distance counts the minimum number of insertions, deletions, and substitutions needed to transform one string into another.

Each match returns:

- **Position** -- start and end coordinates in the target sequence
- **Cost** -- the edit distance (number of mismatches + indels)
- **CIGAR string** -- alignment representation
- **Strand** -- which strand the match was found on

## IUPAC Ambiguity Handling

Sassy supports IUPAC ambiguity codes natively through its `Iupac` profile. Ambiguous bases in primer sequences are expanded during matching -- for example, an `R` in a primer matches both `A` and `G` in the target without counting as a mismatch.

This is handled at the SIMD level, so ambiguous primers do not incur a performance penalty compared to standard ACGT primers.

See [Primer Formats](primer-formats.md#iupac-ambiguity-codes) for the full IUPAC code table.

## Forward/Reverse Primer Pairing

The core PCR algorithm works as follows:

1. **Search for forward primer** -- find all approximate matches of the forward primer in the (possibly extended) sequence

2. **Search for reverse primer downstream** -- for each forward match, search for the reverse complement of the reverse primer in the region from the forward match end to `forward_start + max_len`

3. **Validate product** -- check that the product satisfies:
    - Length between `--min-len` and `--max-len`
    - N-base fraction below `--max-n-fraction`
    - Identity above `--min-identity` for both primers

4. **Extract sequence** -- extract the product from the sequence, either including primers (default) or excluding them (`--trim-primers`)

!!! note "Reverse primer orientation"
    The reverse primer is specified in 5'-to-3' direction (as you would order it). Amplirust internally converts it to its reverse complement to search for the binding site on the forward strand.

## Circular Genome Handling

When `--circular` is enabled, Amplirust extends the sequence to detect products that wrap around the origin of replication.

The extension strategy appends the first `min(max_product_len - 1, sequence_length)` bases of the sequence to its end:

```
Original:  [====SEQUENCE====]
Extended:  [====SEQUENCE====][====]
                              ↑ first N bases repeated
```

This allows the algorithm to find forward primer matches near the end and reverse primer matches that fall in the repeated region. Products that span the junction are flagged with the `_wrap` suffix and `is_circular_wrap = true` in the TSV output.

Positions are mapped back to the original coordinate system using modulo arithmetic.

## Reverse Complement Strand Search

When `--search-rc` is enabled, Amplirust searches both DNA strands:

1. **Forward strand** -- the sequence as provided
2. **Reverse complement strand** -- the full reverse complement of the sequence

The same primer pairing algorithm runs on both strands. Products found on the reverse complement strand are:

- Output in their reverse complement orientation
- Labeled with `strand=-` and the `_rc` suffix in the header
- Position coordinates mapped back to the original forward strand

## Edit Distance vs Identity Filtering

Amplirust provides two complementary filtering mechanisms:

### Edit Distance (`-k, --max-errors`)

Controls the search itself. Only matches within this edit distance are returned by the SIMD search engine. This is the primary performance knob -- lower values search faster.

**Default:** 2

### Identity Threshold (`--min-identity`)

A post-filter applied after matching. Identity is calculated as:

```
identity = (alignment_length - edit_distance) / alignment_length
```

For example, a 20 bp primer with 2 mismatches has identity = 18/20 = 0.90.

**Default:** 0.0 (disabled)

!!! tip "When to use each"
    Use `-k` to control computational cost (lower = faster). Use `--min-identity` for biological relevance when primer lengths vary -- a 2-mismatch hit in a 15 bp primer (87% identity) is less reliable than in a 25 bp primer (92% identity).
