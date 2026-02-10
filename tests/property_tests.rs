use proptest::prelude::*;

use amplirust::utils::{
    calculate_identity, circular_to_original_pos, is_circular_wrap, make_circular_searchable,
};

/// Strategy to generate random DNA sequences of length in `len_range`.
fn dna_sequence(len_range: std::ops::Range<usize>) -> impl Strategy<Value = Vec<u8>> {
    proptest::collection::vec(
        proptest::sample::select(vec![b'A', b'C', b'G', b'T']),
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
}
