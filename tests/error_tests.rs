//! Error path and edge case integration tests for amplirust.
//!
//! TEST-05: Malformed input handling (no panics, graceful errors)
//! TEST-06: Identity filtering thresholds
//! TEST-07: Extreme primer edge cases

use std::path::{Path, PathBuf};

use amplirust::input::{SequenceRecord, read_sequences_from_file};
use amplirust::matcher::MatchConfig;
use amplirust::output::validate_output_writable;
use amplirust::pcr::{PcrConfig, find_pcr_products};
use amplirust::primer::{PrimerPair, parse_primers};
use amplirust::utils::reverse_complement;

fn test_data_path(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(filename)
}

// ============================================================================
// TEST-05: Malformed input handling
// ============================================================================

#[test]
fn test_malformed_fasta_no_panic() {
    // FASTA with no > headers -- just bare sequence lines.
    // The primary assertion is that this line runs without panicking.
    let result = read_sequences_from_file(&test_data_path("malformed.fasta"), 0);

    // Parser should handle gracefully: either Ok with no valid records, or Err.
    if let Ok(records) = result {
        // If Ok, records should be empty or contain only headerless entries.
        // seq_io returns zero records for headerless FASTA (no '>' seen).
        assert!(
            records.is_empty()
                || records
                    .iter()
                    .all(|r| r.header.is_empty() || r.sequence.is_empty()),
            "Malformed FASTA should produce no valid records, got {} record(s)",
            records.len()
        );
    }
    // An Err is also acceptable -- validation caught the problem.
}

#[test]
fn test_truncated_genbank_no_panic() {
    // GenBank file truncated mid-ORIGIN section (no // terminator).
    // Per Phase 2 decision: truncated GenBank emits warning but returns partial data.
    let result = read_sequences_from_file(&test_data_path("truncated.gb"), 0);

    // Must not panic. Ok with partial data is the expected outcome.
    if let Ok(records) = result {
        // Truncated GenBank may return partial sequence data.
        // The parser should handle the missing // gracefully.
        if !records.is_empty() {
            // If records returned, verify they contain some sequence data.
            assert!(
                records.iter().any(|r| !r.sequence.is_empty()),
                "Truncated GenBank records should contain partial sequence data"
            );
        }
        // Empty records is also acceptable if parser skips incomplete records.
    }
    // An Err is also acceptable for severely truncated files.
}

#[test]
fn test_bad_primers_csv_returns_error() {
    // CSV with only 2 columns (name, forward) instead of required 3 (name, forward, reverse).
    let path = test_data_path("bad_primers.csv");
    let result = parse_primers(&path.to_string_lossy());

    assert!(
        result.is_err(),
        "Bad CSV with wrong column count should return error"
    );

    // Error message should contain diagnostic information about the column count.
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("column") || err_msg.contains('2') || err_msg.contains('3'),
        "Error should mention column count issue, got: {err_msg}"
    );
}

#[test]
fn test_unwritable_output_returns_error() {
    // Path in a nonexistent directory.
    let result = validate_output_writable(Path::new("/proc/nonexistent/dir/output.fasta"));

    assert!(
        result.is_err(),
        "Unwritable output path should return error"
    );

    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("does not exist") || err_msg.contains("not writable"),
        "Error should indicate directory problem, got: {err_msg}"
    );
}

#[test]
fn test_empty_fasta_returns_empty_records() {
    // Create a temp file with just whitespace, suffix ".fasta"
    let temp = tempfile::Builder::new()
        .suffix(".fasta")
        .tempfile()
        .expect("Failed to create temp file");
    std::fs::write(temp.path(), "   \n\n  \n").expect("Failed to write temp file");

    let result = read_sequences_from_file(temp.path(), 0);

    // Must not panic. Ok with empty vec is expected.
    if let Ok(records) = result {
        assert!(
            records.is_empty(),
            "Whitespace-only FASTA should produce no records, got {}",
            records.len()
        );
    }
    // An Err is also acceptable for files with no valid content.
}

#[test]
fn test_empty_sequence_no_panic() {
    // SequenceRecord with empty sequence vec.
    let record = SequenceRecord {
        header: "empty".to_string(),
        sequence: Vec::new(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer =
        PrimerPair::new("test", b"ACGTACGT", b"TGCATGCA").expect("Valid primer should construct");

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 4,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    // Must not panic. Should return empty products.
    let products = find_pcr_products(&record, &primer, &config);
    assert!(
        products.is_empty(),
        "Empty sequence should produce no products"
    );
}

// ============================================================================
// TEST-06: Identity filtering
// ============================================================================

/// Build a 200bp sequence with a mutated forward primer at position 50 and
/// exact reverse primer RC at position 150.
///
/// Forward primer: ACGTACGTACGTACGTACGT (20bp)
/// Mutated fwd at pos 50: ACGTAAGTAAGTAAGTACGT (4 substitutions => 80% identity)
/// Reverse primer: TGCATGCATGCATGCATGCA (20bp), RC = TGCATGCATGCATGCATGCA
/// (TGCA reversed = ACGT, complemented each: T->A, G->C, C->G, A->T => ACGT...
///  Actually: reverse of TGCATGCATGCATGCATGCA = ACGTACGTACGTACGTACGT,
///  complement = TGCATGCATGCATGCATGCA. So RC = TGCATGCATGCATGCATGCA -- palindromic!)
///
/// To avoid palindromic issues, use asymmetric primers:
/// Forward: AACGAACGAACGAACGAACG (20bp)
/// Reverse: TTGCTTGCTTGCTTGCTTGC (20bp), RC = GCAAGCAAGCAAGCAAGCAA
fn build_identity_test_sequence() -> (Vec<u8>, PrimerPair) {
    let fwd_primer = b"AACGAACGAACGAACGAACG"; // 20bp
    let rev_primer = b"TTGCTTGCTTGCTTGCTTGC"; // 20bp
    let rev_primer_rc = reverse_complement(rev_primer); // What we embed in the sequence

    // Mutated forward: 4 substitutions at positions 4,5,9,10 (0-indexed)
    // Original: AACGAACGAACGAACGAACG
    // Mutated:  AACGTTCGTTCGAACGAACG (4 changes: A->T at 4,5 and A->T at 8,9)
    // Actually let me be more careful:
    //   pos: 0123456789...
    //   orig: A A C G A A C G A A C G A A C G A A C G
    //   mut:  A A C G T T C G T T C G A A C G A A C G
    // Changes at positions 4,5,8,9: A->T (4 substitutions = 80% identity for 20bp)
    let mutated_fwd = b"AACGTTCGTTCGAACGAACG";

    // Build 200bp sequence
    let mut seq = Vec::with_capacity(200);
    // Padding: 50 bp of A's
    seq.extend_from_slice(&[b'A'; 50]);
    // Mutated forward primer at position 50
    seq.extend_from_slice(mutated_fwd);
    // Filler: 80 bp of T's (positions 70-149)
    seq.extend_from_slice(&[b'T'; 80]);
    // Exact reverse primer RC at position 150
    seq.extend_from_slice(&rev_primer_rc);
    // Padding: remaining to fill 200bp (200 - 50 - 20 - 80 - 20 = 30)
    seq.extend_from_slice(&[b'A'; 30]);

    assert_eq!(seq.len(), 200, "Test sequence should be 200bp");

    let primer =
        PrimerPair::new("identity_test", fwd_primer, rev_primer).expect("Valid primer pair");
    (seq, primer)
}

#[test]
fn test_identity_filter_rejects_low_identity_match() {
    let (seq, primer) = build_identity_test_sequence();

    let record = SequenceRecord {
        header: "identity_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("test.fasta"),
    };

    // 4 mismatches in 20bp = 80% identity.
    // max_errors=4 allows the alignment to be found by sassy,
    // but min_identity=0.9 (90% threshold) should filter it out.
    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 4,
            min_identity: 0.9,
            search_rc: false,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    // 80% identity forward match should be filtered by 90% threshold
    assert!(
        products.is_empty(),
        "4 mismatches in 20bp = 80% identity, should be filtered by 90% threshold. Got {} product(s)",
        products.len()
    );
}

#[test]
fn test_identity_filter_passes_same_match_with_lower_threshold() {
    let (seq, primer) = build_identity_test_sequence();

    let record = SequenceRecord {
        header: "identity_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("test.fasta"),
    };

    // Same sequence, same primers, but min_identity=0.7 (70% threshold).
    // 80% identity should now pass.
    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 4,
            min_identity: 0.7,
            search_rc: false,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    // 80% identity passes 70% threshold -- positive control
    assert!(
        !products.is_empty(),
        "4 mismatches in 20bp = 80% identity, should pass 70% threshold"
    );
}

#[test]
fn test_identity_filter_exact_match_always_passes() {
    // Sanity check: exact match with strictest threshold should always pass.
    let fwd = b"ACGTACGTACGTACGTACGT"; // 20bp
    let rev = b"TGCATGCATGCATGCATGCA"; // 20bp
    let rev_rc = reverse_complement(rev);

    // Build sequence with exact forward and exact reverse RC
    let mut seq = Vec::with_capacity(120);
    seq.extend_from_slice(&[b'A'; 20]);
    seq.extend_from_slice(fwd); // 20bp forward at pos 20
    seq.extend_from_slice(&[b'T'; 40]); // filler
    seq.extend_from_slice(&rev_rc); // 20bp reverse RC at pos 80
    // Pad to at least 100bp
    while seq.len() < 120 {
        seq.push(b'A');
    }

    let record = SequenceRecord {
        header: "exact_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("exact", fwd, rev).expect("Valid primer pair");

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0, // 100% identity required
            search_rc: false,
        },
        min_len: 20,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    assert!(
        !products.is_empty(),
        "Exact match with 100% identity threshold should find product"
    );
}

// ============================================================================
// TEST-07: Edge-case primers
// ============================================================================

#[test]
fn test_very_short_primer_3bp() {
    // 3bp primers are below the recommended minimum (10bp) but should not panic.
    let primer_result = PrimerPair::new("short", b"ACG", b"TGC");

    if let Ok(primer) = primer_result {
        // RC of TGC = GCA
        let mut seq = Vec::with_capacity(50);
        seq.extend_from_slice(b"AAAAA");
        seq.extend_from_slice(b"ACG"); // forward primer
        seq.extend_from_slice(b"NNNNNNNNNN");
        seq.extend_from_slice(&reverse_complement(&primer.reverse)); // RC of reverse
        seq.extend_from_slice(b"AAAAA");

        let record = SequenceRecord {
            header: "short_primer_test".to_string(),
            sequence: seq,
            source_file: PathBuf::from("test.fasta"),
        };

        let config = PcrConfig {
            match_config: MatchConfig {
                max_errors: 0,
                min_identity: 1.0,
                search_rc: false,
            },
            min_len: 3,
            max_len: 100,
            circular: false,
            trim_primers: false,
            max_n_fraction: 1.0,
        };

        // Key assertion: no panic. Products may or may not be found.
        let _products = find_pcr_products(&record, &primer, &config);
    }
    // If PrimerPair rejects 3bp primers, that is also acceptable.
}

#[test]
fn test_very_long_primer_60bp() {
    // Generate a deterministic 60bp primer from repeating pattern.
    let fwd_60: Vec<u8> = b"ACGTACGTACGTACGTACGT"
        .iter()
        .cycle()
        .take(60)
        .copied()
        .collect();
    let rev_60: Vec<u8> = b"TGCATGCATGCATGCATGCA"
        .iter()
        .cycle()
        .take(60)
        .copied()
        .collect();
    let rev_60_rc = reverse_complement(&rev_60);

    let primer = PrimerPair::new("long60", &fwd_60, &rev_60).expect("60bp primer should be valid");

    // Build sequence with the 60bp primers embedded
    let mut seq = Vec::with_capacity(300);
    seq.extend_from_slice(&[b'A'; 20]);
    seq.extend_from_slice(&fwd_60); // 60bp forward primer at pos 20
    seq.extend_from_slice(&[b'T'; 100]); // 100bp filler
    seq.extend_from_slice(&rev_60_rc); // 60bp reverse RC at pos 180
    seq.extend_from_slice(&[b'A'; 20]);

    let record = SequenceRecord {
        header: "long_primer_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("test.fasta"),
    };

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 100,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    // Key assertion: no panic. If product found, verify coordinates are sensible.
    let products = find_pcr_products(&record, &primer, &config);
    for product in &products {
        assert!(
            product.full_length >= 100,
            "Product with 60bp primers should be at least 100bp"
        );
    }
}

#[test]
fn test_all_n_primer() {
    // 10bp of all N -- N is valid IUPAC (matches any base).
    let primer_result = PrimerPair::new("all_n", b"NNNNNNNNNN", b"NNNNNNNNNN");

    if let Ok(primer) = primer_result {
        let record = SequenceRecord {
            header: "all_n_test".to_string(),
            sequence: b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), // 32bp
            source_file: PathBuf::from("test.fasta"),
        };

        let config = PcrConfig {
            match_config: MatchConfig {
                max_errors: 0,
                min_identity: 0.0, // Very permissive -- N matches anything in IUPAC
                search_rc: false,
            },
            min_len: 10,
            max_len: 100,
            circular: false,
            trim_primers: false,
            max_n_fraction: 1.0,
        };

        // Key assertion: no panic. May find many products or none depending
        // on how sassy handles all-N patterns.
        let _products = find_pcr_products(&record, &primer, &config);
    }
    // If PrimerPair rejects all-N primers, that is also acceptable.
}

#[test]
fn test_palindromic_primer() {
    // ACGTACGT is a palindromic sequence (RC of ACGTACGT = ACGTACGT).
    // Using the same sequence for both forward and reverse creates a
    // self-complementary primer pair.
    let primer = PrimerPair::new("palindrome", b"ACGTACGT", b"ACGTACGT")
        .expect("Palindromic primer should be valid");

    // Embed ACGTACGT at two positions for forward and reverse matches.
    // RC of ACGTACGT = ACGTACGT, so the matcher searches for ACGTACGT downstream.
    let mut seq = Vec::with_capacity(60);
    seq.extend_from_slice(b"AAAA");
    seq.extend_from_slice(b"ACGTACGT"); // position 4: first occurrence
    seq.extend_from_slice(b"TTTTTTTTTTTT"); // 12bp filler
    seq.extend_from_slice(b"ACGTACGT"); // position 24: second occurrence
    seq.extend_from_slice(b"AAAA");

    let record = SequenceRecord {
        header: "palindrome_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("test.fasta"),
    };

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 4,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    // Key assertion: no panic. Palindromic primers should not cause issues.
    let products = find_pcr_products(&record, &primer, &config);

    // Should find at least one product: from first ACGTACGT to second ACGTACGT.
    assert!(
        !products.is_empty(),
        "Palindromic primer should find at least one product"
    );
}
