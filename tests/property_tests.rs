use proptest::prelude::*;

use amplirust::utils::{
    calculate_identity, circular_to_original_pos, complement, is_circular_wrap,
    make_circular_searchable, reverse_complement,
};

/// Strategy to generate random DNA sequences of length in `len_range`.
fn dna_sequence(len_range: std::ops::Range<usize>) -> impl Strategy<Value = Vec<u8>> {
    proptest::collection::vec(
        proptest::sample::select(vec![b'A', b'C', b'G', b'T']),
        len_range,
    )
}

/// Strategy to generate random IUPAC sequences (all 15 valid ambiguity codes).
fn iupac_sequence(len_range: std::ops::Range<usize>) -> impl Strategy<Value = Vec<u8>> {
    proptest::collection::vec(
        proptest::sample::select(vec![
            b'A', b'C', b'G', b'T', b'R', b'Y', b'S', b'W', b'K', b'M', b'B', b'D', b'H', b'V',
            b'N',
        ]),
        len_range,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    // ── TEST-08: Coordinate math invariants ──────────────────────────────

    /// circular_to_original_pos always returns a value strictly less than original_len.
    #[test]
    fn circular_pos_within_bounds(pos in 0..100_000_usize, original_len in 1..10_000_usize) {
        let result = circular_to_original_pos(pos, original_len);
        prop_assert!(
            result < original_len,
            "circular_to_original_pos({}, {}) = {} which is >= original_len",
            pos,
            original_len,
            result
        );
    }

    /// Non-wrapping positions map to themselves.
    #[test]
    fn circular_wrap_consistency(pos in 0..100_000_usize, original_len in 1..10_000_usize) {
        if !is_circular_wrap(pos, original_len) {
            let result = circular_to_original_pos(pos, original_len);
            prop_assert_eq!(
                result,
                pos,
                "Non-wrapping pos {} should map to itself for original_len {}",
                pos,
                original_len
            );
        }
    }

    /// make_circular_searchable extends by at most max_product_len.saturating_sub(1).min(seq.len()).
    #[test]
    fn circular_extension_bounded(
        seq in dna_sequence(1..500),
        max_product_len in 0..1000_usize
    ) {
        let extended = make_circular_searchable(&seq, max_product_len);
        let extension = extended.len() - seq.len();
        let max_extension = max_product_len.saturating_sub(1).min(seq.len());
        prop_assert!(
            extension <= max_extension,
            "Extension {} exceeds bound {} for seq len {} and max_product_len {}",
            extension,
            max_extension,
            seq.len(),
            max_product_len
        );
    }

    /// calculate_identity always returns a value in [0.0, 1.0].
    #[test]
    fn identity_bounded_zero_one(
        edit_distance in 0..1000_usize,
        alignment_len in 0..1000_usize
    ) {
        let identity = calculate_identity(edit_distance, alignment_len);
        prop_assert!(
            (0.0..=1.0).contains(&identity),
            "calculate_identity({}, {}) = {} which is outside [0.0, 1.0]",
            edit_distance,
            alignment_len,
            identity
        );
    }

    // ── TEST-09: Reverse complement invariants ───────────────────────────

    /// Applying reverse complement twice to any DNA sequence returns the original.
    #[test]
    fn rc_involution_dna(seq in dna_sequence(0..500)) {
        let len = seq.len();
        let rc = reverse_complement(&seq);
        let rc_rc = reverse_complement(&rc);
        prop_assert_eq!(
            rc_rc,
            seq,
            "RC involution failed for DNA sequence of length {}",
            len
        );
    }

    /// RC involution holds for all 15 valid IUPAC ambiguity codes.
    #[test]
    fn rc_involution_iupac(seq in iupac_sequence(0..500)) {
        let len = seq.len();
        let rc = reverse_complement(&seq);
        let rc_rc = reverse_complement(&rc);
        prop_assert_eq!(
            rc_rc,
            seq,
            "RC involution failed for IUPAC sequence of length {}",
            len
        );
    }

    /// Reverse complement always preserves sequence length.
    #[test]
    fn rc_preserves_length(seq in dna_sequence(0..1000)) {
        let rc = reverse_complement(&seq);
        prop_assert_eq!(
            rc.len(),
            seq.len(),
            "RC changed length: {} -> {}",
            seq.len(),
            rc.len()
        );
    }

    /// Single-base complement is its own inverse for all valid IUPAC bases.
    #[test]
    fn rc_complement_involution_single_base(
        base in proptest::sample::select(vec![
            b'A', b'C', b'G', b'T', b'R', b'Y', b'S', b'W', b'K', b'M',
            b'B', b'D', b'H', b'V', b'N',
        ])
    ) {
        let result = complement(complement(base));
        prop_assert_eq!(
            result,
            base,
            "complement(complement({})) = {} (expected {})",
            base as char,
            result as char,
            base as char
        );
    }
}
