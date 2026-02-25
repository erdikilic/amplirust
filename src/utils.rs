//! DNA sequence utilities: complement, reverse complement, and circular genome helpers.
//!
//! Provides IUPAC-aware nucleotide complementation, circular sequence extension
//! for wrap-around product detection, and coordinate mapping functions.

/// DNA complement table (handles IUPAC codes)
/// Index by ASCII value to get complement
const COMPLEMENT_TABLE: [u8; 128] = {
    let mut table = [b'N'; 128];

    // Standard bases
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'C' as usize] = b'G';
    table[b'U' as usize] = b'A'; // RNA

    // Lowercase
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'g' as usize] = b'c';
    table[b'c' as usize] = b'g';
    table[b'u' as usize] = b'a';

    // IUPAC ambiguity codes
    table[b'R' as usize] = b'Y'; // R = A|G -> Y = T|C
    table[b'Y' as usize] = b'R'; // Y = T|C -> R = A|G
    table[b'S' as usize] = b'S'; // S = G|C -> S = C|G (self-complementary)
    table[b'W' as usize] = b'W'; // W = A|T -> W = T|A (self-complementary)
    table[b'K' as usize] = b'M'; // K = G|T -> M = C|A
    table[b'M' as usize] = b'K'; // M = A|C -> K = T|G
    table[b'B' as usize] = b'V'; // B = C|G|T -> V = G|C|A
    table[b'V' as usize] = b'B'; // V = A|C|G -> B = T|G|C
    table[b'D' as usize] = b'H'; // D = A|G|T -> H = T|C|A
    table[b'H' as usize] = b'D'; // H = A|C|T -> D = T|G|A
    table[b'N' as usize] = b'N'; // N = any -> N

    // Lowercase IUPAC
    table[b'r' as usize] = b'y';
    table[b'y' as usize] = b'r';
    table[b's' as usize] = b's';
    table[b'w' as usize] = b'w';
    table[b'k' as usize] = b'm';
    table[b'm' as usize] = b'k';
    table[b'b' as usize] = b'v';
    table[b'v' as usize] = b'b';
    table[b'd' as usize] = b'h';
    table[b'h' as usize] = b'd';
    table[b'n' as usize] = b'n';

    table
};

/// Get the complement of a single nucleotide
#[inline]
#[must_use]
pub fn complement(base: u8) -> u8 {
    if base < 128 {
        COMPLEMENT_TABLE[base as usize]
    } else {
        b'N'
    }
}

/// Compute the reverse complement of a DNA sequence
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Compute reverse complement in place (for efficiency)
pub fn reverse_complement_into(seq: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.reserve(seq.len());
    for &b in seq.iter().rev() {
        output.push(complement(b));
    }
}

/// Calculate percentage identity from edit distance and alignment length
/// Returns value between 0.0 and 1.0
#[must_use]
pub fn calculate_identity(edit_distance: usize, alignment_len: usize) -> f64 {
    if alignment_len == 0 {
        return 0.0;
    }
    let matches = alignment_len.saturating_sub(edit_distance);
    matches as f64 / alignment_len as f64
}

/// Extend a sequence for circular genome searching
/// Appends the beginning of the sequence to allow wrap-around matches
#[must_use]
pub fn make_circular_searchable(seq: &[u8], max_product_len: usize) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }

    // We need to extend by at most max_product_len - 1 to catch wrap-around
    // But also don't extend more than the sequence length itself
    let extend_len = (max_product_len.saturating_sub(1)).min(seq.len());

    let mut extended = Vec::with_capacity(seq.len() + extend_len);
    extended.extend_from_slice(seq);
    extended.extend_from_slice(&seq[..extend_len]);
    extended
}

/// Check if a position in an extended circular sequence wraps around
#[must_use]
pub fn is_circular_wrap(pos: usize, original_len: usize) -> bool {
    pos >= original_len
}

/// Convert a position from extended circular sequence back to original coordinates
#[must_use]
pub fn circular_to_original_pos(pos: usize, original_len: usize) -> usize {
    pos % original_len
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── complement ──────────────────────────────────────────────────────

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_complement_iupac() {
        assert_eq!(complement(b'R'), b'Y');
        assert_eq!(complement(b'Y'), b'R');
        assert_eq!(complement(b'S'), b'S');
        assert_eq!(complement(b'W'), b'W');
        assert_eq!(complement(b'K'), b'M');
        assert_eq!(complement(b'M'), b'K');
    }

    #[test]
    fn test_complement_rna() {
        assert_eq!(complement(b'U'), b'A');
        assert_eq!(complement(b'u'), b'a');
    }

    #[test]
    fn test_complement_non_ascii() {
        // Bytes >= 128 should return N
        assert_eq!(complement(128), b'N');
        assert_eq!(complement(255), b'N');
        assert_eq!(complement(200), b'N');
    }

    #[test]
    fn test_complement_lowercase() {
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(b't'), b'a');
        assert_eq!(complement(b'g'), b'c');
        assert_eq!(complement(b'c'), b'g');
    }

    // ── reverse_complement ──────────────────────────────────────────────

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_rc_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_rc_single_base() {
        assert_eq!(reverse_complement(b"A"), b"T");
    }

    #[test]
    fn test_rc_palindrome() {
        // ACGT is its own reverse complement
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
    }

    #[test]
    fn test_rc_all_same() {
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }

    #[test]
    fn test_rc_iupac_codes() {
        // R->Y, Y->R, S->S, W->W, K->M, M->K
        // Reverse of RYSWKM = MKWSYR, then complement each:
        // M->K, K->M, W->W, S->S, Y->R, R->Y => KMWSRY
        assert_eq!(reverse_complement(b"RYSWKM"), b"KMWSRY");
    }

    #[test]
    fn test_rc_lowercase() {
        assert_eq!(reverse_complement(b"acgt"), b"acgt");
    }

    #[test]
    fn test_rc_involution() {
        // Applying reverse complement twice should return the original
        let sequences: &[&[u8]] = &[b"ACGT", b"RYSWKM", b"ATCGATCG", b"A", b""];
        for seq in sequences {
            let rc = reverse_complement(seq);
            let rc_rc = reverse_complement(&rc);
            assert_eq!(
                &rc_rc,
                *seq,
                "RC involution failed for {:?}",
                String::from_utf8_lossy(seq)
            );
        }
    }

    // ── calculate_identity ──────────────────────────────────────────────

    #[test]
    fn test_calculate_identity() {
        assert!((calculate_identity(0, 20) - 1.0).abs() < 0.001);
        assert!((calculate_identity(2, 20) - 0.9).abs() < 0.001);
        assert!((calculate_identity(10, 20) - 0.5).abs() < 0.001);
        assert!((calculate_identity(0, 0) - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_identity_perfect() {
        assert!((calculate_identity(0, 20) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_identity_zero_len() {
        assert!((calculate_identity(0, 0) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_identity_half() {
        assert!((calculate_identity(10, 20) - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_identity_distance_exceeds_len() {
        // saturating_sub prevents negative: (20 - 25).saturating = 0
        assert!((calculate_identity(25, 20) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_identity_distance_equals_len() {
        assert!((calculate_identity(20, 20) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_identity_one_error() {
        assert!((calculate_identity(1, 100) - 0.99).abs() < f64::EPSILON);
    }

    // ── make_circular_searchable ────────────────────────────────────────

    #[test]
    fn test_circular_searchable() {
        let seq = b"ABCDEFGH";
        let extended = make_circular_searchable(seq, 4);
        assert_eq!(&extended, b"ABCDEFGHABC");
    }

    #[test]
    fn test_make_circular_searchable_empty() {
        assert_eq!(make_circular_searchable(b"", 10), b"");
    }

    #[test]
    fn test_make_circular_searchable_max_len_zero() {
        // max_product_len=0 => extend_len = 0.saturating_sub(1) = 0
        let seq = b"ACGT";
        let result = make_circular_searchable(seq, 0);
        assert_eq!(&result, b"ACGT");
    }

    #[test]
    fn test_make_circular_searchable_max_len_one() {
        // max_product_len=1 => extend_len = (1-1).min(4) = 0
        let seq = b"ACGT";
        let result = make_circular_searchable(seq, 1);
        assert_eq!(&result, b"ACGT");
    }

    #[test]
    fn test_make_circular_searchable_seq_shorter_than_max() {
        // seq len 4, max_product_len=10 => extend_len = (10-1).min(4) = 4
        let seq = b"ABCD";
        let result = make_circular_searchable(seq, 10);
        assert_eq!(&result, b"ABCDABCD");
    }

    #[test]
    fn test_make_circular_searchable_seq_equal_max() {
        // seq len 8, max_product_len=8 => extend_len = (8-1).min(8) = 7
        let seq = b"ABCDEFGH";
        let result = make_circular_searchable(seq, 8);
        assert_eq!(&result, b"ABCDEFGHABCDEFG");
    }

    #[test]
    fn test_make_circular_searchable_seq_longer_than_max() {
        // seq len 20, max_product_len=4 => extend_len = (4-1).min(20) = 3
        let seq = b"ABCDEFGHIJKLMNOPQRST";
        let result = make_circular_searchable(seq, 4);
        assert_eq!(&result, b"ABCDEFGHIJKLMNOPQRSTABC");
    }

    #[test]
    fn test_make_circular_searchable_single_base() {
        // seq=[A], max_product_len=5 => extend_len = (5-1).min(1) = 1
        let seq = b"A";
        let result = make_circular_searchable(seq, 5);
        assert_eq!(&result, b"AA");
    }

    // ── is_circular_wrap ────────────────────────────────────────────────

    #[test]
    fn test_circular_wrap() {
        assert!(!is_circular_wrap(5, 10));
        assert!(is_circular_wrap(10, 10));
        assert!(is_circular_wrap(15, 10));
    }

    #[test]
    fn test_circular_wrap_at_zero() {
        assert!(!is_circular_wrap(0, 10));
    }

    #[test]
    fn test_circular_wrap_at_last() {
        assert!(!is_circular_wrap(9, 10));
    }

    #[test]
    fn test_circular_wrap_at_boundary() {
        assert!(is_circular_wrap(10, 10));
    }

    #[test]
    fn test_circular_wrap_well_beyond() {
        assert!(is_circular_wrap(25, 10));
    }

    #[test]
    fn test_circular_wrap_original_len_one() {
        assert!(!is_circular_wrap(0, 1));
        assert!(is_circular_wrap(1, 1));
    }

    // ── circular_to_original_pos ────────────────────────────────────────

    #[test]
    fn test_circular_to_original() {
        assert_eq!(circular_to_original_pos(5, 10), 5);
        assert_eq!(circular_to_original_pos(10, 10), 0);
        assert_eq!(circular_to_original_pos(12, 10), 2);
    }

    #[test]
    fn test_circular_pos_zero() {
        assert_eq!(circular_to_original_pos(0, 10), 0);
        assert_eq!(circular_to_original_pos(0, 1), 0);
        assert_eq!(circular_to_original_pos(0, 100), 0);
    }

    #[test]
    fn test_circular_pos_at_len() {
        // pos=10, original_len=10 => 10 % 10 = 0 (wraps)
        assert_eq!(circular_to_original_pos(10, 10), 0);
    }

    #[test]
    fn test_circular_pos_beyond_len() {
        // pos=13, original_len=10 => 13 % 10 = 3
        assert_eq!(circular_to_original_pos(13, 10), 3);
    }

    #[test]
    fn test_circular_pos_original_len_one() {
        // Everything maps to 0 when original_len=1
        assert_eq!(circular_to_original_pos(0, 1), 0);
        assert_eq!(circular_to_original_pos(1, 1), 0);
        assert_eq!(circular_to_original_pos(5, 1), 0);
    }

    #[test]
    fn test_circular_pos_double_wrap() {
        // pos=25, original_len=10 => 25 % 10 = 5
        assert_eq!(circular_to_original_pos(25, 10), 5);
    }

    // ── reverse_complement_into (mutation gap) ────────────────────────────

    #[test]
    fn test_reverse_complement_into_basic() {
        let mut output = Vec::new();
        reverse_complement_into(b"ACGT", &mut output);
        assert_eq!(output, b"ACGT"); // ACGT is palindromic RC
    }

    #[test]
    fn test_reverse_complement_into_non_palindrome() {
        let mut output = Vec::new();
        reverse_complement_into(b"AAAA", &mut output);
        assert_eq!(output, b"TTTT");
    }

    #[test]
    fn test_reverse_complement_into_agrees_with_reverse_complement() {
        let sequences: &[&[u8]] = &[b"ATCG", b"RYSWKM", b"A", b"GATTACA"];
        for seq in sequences {
            let mut output = Vec::new();
            reverse_complement_into(seq, &mut output);
            assert_eq!(
                output,
                reverse_complement(seq),
                "reverse_complement_into disagrees for {:?}",
                String::from_utf8_lossy(seq)
            );
        }
    }

    #[test]
    fn test_reverse_complement_into_empty() {
        let mut output = Vec::new();
        reverse_complement_into(b"", &mut output);
        assert!(output.is_empty());
    }

    #[test]
    fn test_reverse_complement_into_clears_output() {
        let mut output = vec![b'X', b'Y', b'Z'];
        reverse_complement_into(b"AT", &mut output);
        assert_eq!(output, b"AT"); // AT -> RC = AT; old contents cleared
    }
}
