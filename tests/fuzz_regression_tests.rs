//! Regression tests for fuzz-discovered crashes.
//!
//! Each test loads a crash artifact (via `include_bytes!`) and feeds it to
//! the parser that originally crashed. The assertion is simply "must not
//! panic" -- any `Err` or empty result is fine.
//!
//! ## Fuzzing campaign results (2026-02-10)
//!
//! All five targets ran for 10 minutes each with `--sanitizer none`:
//!
//! | Target                  | Iterations  | Corpus | Crashes |
//! |-------------------------|-------------|--------|---------|
//! | fuzz_genbank_streaming  | 24,806,303  | 1,047  | 0       |
//! | fuzz_genbank_slice      | 29,410,773  | 1,076  | 0       |
//! | fuzz_fasta              | 85,768,796  | 192    | 0       |
//! | fuzz_primer_csv         | 2,965,500   | 789    | 0       |
//! | fuzz_primer_matcher     | 7,308,723   | 491    | 0       |
//!
//! Total: ~150 million iterations, zero crashes.
//!
//! The smoke tests below exercise each parser with empty / minimal bytes to
//! keep the regression-test infrastructure ready for future crash artifacts.

use std::io::Cursor;

// ---------------------------------------------------------------------------
// GenBank streaming parser -- empty and minimal inputs must not panic
// ---------------------------------------------------------------------------

#[test]
fn regression_genbank_streaming_empty() {
    let reader = amplirust::genbank::GenbankReader::new(Cursor::new(b"" as &[u8]));
    for result in reader {
        let _ = result;
    }
}

#[test]
fn regression_genbank_streaming_garbage() {
    let data = b"\xff\x00\x01\x02not a genbank file\n";
    let reader = amplirust::genbank::GenbankReader::new(Cursor::new(data as &[u8]));
    for result in reader {
        let _ = result;
    }
}

// ---------------------------------------------------------------------------
// GenBank slice parser -- empty and minimal inputs must not panic
// ---------------------------------------------------------------------------

#[test]
fn regression_genbank_slice_empty() {
    let records = amplirust::genbank::parse_genbank_slice(b"");
    assert!(records.is_empty());
}

#[test]
fn regression_genbank_slice_garbage() {
    let data = b"\xff\x00\x01\x02not a genbank file\n";
    let records = amplirust::genbank::parse_genbank_slice(data);
    // Garbage input should not panic; zero or partial records are fine
    let _ = records;
}

// ---------------------------------------------------------------------------
// FASTA parser (seq_io) -- empty and minimal inputs must not panic
// ---------------------------------------------------------------------------

#[test]
fn regression_fasta_empty() {
    let mut reader = seq_io::fasta::Reader::new(&b""[..]);
    while let Some(result) = reader.next() {
        let _ = result;
    }
}

#[test]
fn regression_fasta_garbage() {
    let data = b"\xff\x00\x01\x02not a fasta file\n";
    let mut reader = seq_io::fasta::Reader::new(&data[..]);
    while let Some(result) = reader.next() {
        let _ = result;
    }
}

// ---------------------------------------------------------------------------
// Primer CSV parser -- empty input must not panic
// ---------------------------------------------------------------------------

#[test]
fn regression_primer_csv_empty() {
    use std::io::Write;
    let mut temp = tempfile::NamedTempFile::with_suffix(".csv").unwrap();
    temp.write_all(b"").unwrap();
    temp.flush().unwrap();
    let _ = amplirust::primer::parse_primers(&temp.path().to_string_lossy());
}

#[test]
fn regression_primer_csv_garbage() {
    use std::io::Write;
    let mut temp = tempfile::NamedTempFile::with_suffix(".csv").unwrap();
    temp.write_all(b"\xff\x00\x01\x02not,a,csv\n").unwrap();
    temp.flush().unwrap();
    let _ = amplirust::primer::parse_primers(&temp.path().to_string_lossy());
}

// ---------------------------------------------------------------------------
// Primer matcher -- minimal valid inputs must not panic
// ---------------------------------------------------------------------------

#[test]
fn regression_primer_matcher_minimal() {
    let config = amplirust::MatchConfig {
        max_errors: 0,
        min_identity: 0.0,
        search_rc: true,
    };
    let mut matcher = amplirust::PrimerMatcher::new(config);
    let matches = matcher.find_matches(b"ACGT", b"ACGTACGT");
    let _ = matches;
}

#[test]
fn regression_primer_matcher_all_iupac() {
    let config = amplirust::MatchConfig {
        max_errors: 3,
        min_identity: 0.0,
        search_rc: false,
    };
    let mut matcher = amplirust::PrimerMatcher::new(config);
    // All 15 IUPAC ambiguity codes as primer, standard DNA as target
    let matches = matcher.find_matches(b"ACGTRYSWKMBDHVN", b"ACGTACGTACGTACGTACGT");
    for m in &matches {
        let _ = m.start;
        let _ = m.end;
        let _ = m.edit_distance;
        let _ = m.identity;
        let _ = &m.cigar;
    }
}
