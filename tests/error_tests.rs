//! Error path and edge case integration tests for amplirust.
//!
//! TEST-05: Malformed input handling (no panics, graceful errors)
//! TEST-06: Identity filtering thresholds
//! TEST-07: Extreme primer edge cases

use std::path::{Path, PathBuf};

use amplirust::input::{read_sequences_from_file, SequenceRecord};
use amplirust::matcher::MatchConfig;
use amplirust::output::validate_output_writable;
use amplirust::pcr::{find_pcr_products, PcrConfig};
use amplirust::primer::{parse_primers, PrimerPair};
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
    match result {
        Ok(records) => {
            // If Ok, records should be empty or contain only headerless entries.
            // seq_io returns zero records for headerless FASTA (no '>' seen).
            assert!(
                records.is_empty()
                    || records.iter().all(|r| r.header.is_empty() || r.sequence.is_empty()),
                "Malformed FASTA should produce no valid records, got {} record(s)",
                records.len()
            );
        }
        Err(_) => {
            // An error is also acceptable -- validation caught the problem.
        }
    }
}

#[test]
fn test_truncated_genbank_no_panic() {
    // GenBank file truncated mid-ORIGIN section (no // terminator).
    // Per Phase 2 decision: truncated GenBank emits warning but returns partial data.
    let result = read_sequences_from_file(&test_data_path("truncated.gb"), 0);

    // Must not panic. Ok with partial data is the expected outcome.
    match result {
        Ok(records) => {
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
        Err(_) => {
            // An error is also acceptable for severely truncated files.
        }
    }
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
        err_msg.contains("column") || err_msg.contains("2") || err_msg.contains("3"),
        "Error should mention column count issue, got: {err_msg}"
    );
}

#[test]
fn test_unwritable_output_returns_error() {
    // Path in a nonexistent directory.
    let result = validate_output_writable(Path::new("/proc/nonexistent/dir/output.fasta"));

    assert!(result.is_err(), "Unwritable output path should return error");

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
    match result {
        Ok(records) => {
            assert!(
                records.is_empty(),
                "Whitespace-only FASTA should produce no records, got {}",
                records.len()
            );
        }
        Err(_) => {
            // An error is also acceptable for files with no valid content.
        }
    }
}

#[test]
fn test_empty_sequence_no_panic() {
    // SequenceRecord with empty sequence vec.
    let record = SequenceRecord {
        header: "empty".to_string(),
        sequence: Vec::new(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGTACGT", b"TGCATGCA")
        .expect("Valid primer should construct");

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
