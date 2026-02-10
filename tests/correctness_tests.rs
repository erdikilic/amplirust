//! Biological correctness integration tests for amplirust.
//!
//! These tests prove the core biological correctness claim: every PCR product
//! reported is biologically correct with accurate coordinates and sequence content.
//!
//! TEST-01: 16S primer amplicon at expected coordinates
//! TEST-02: Circular plasmid wrap-around detection
//! TEST-03: Multi-record GenBank independent per-record results
//! TEST-04: RC strand product validation

use std::io::Write;
use std::path::PathBuf;

use tempfile::NamedTempFile;

use amplirust::input::{read_sequences_from_file, SequenceRecord};
use amplirust::matcher::MatchConfig;
use amplirust::pcr::{find_pcr_products, PcrConfig};
use amplirust::primer::{parse_primers, PrimerPair};
use amplirust::utils::reverse_complement;
use sassy::Strand;

// ============================================================================
// Helper functions
// ============================================================================

/// Return the path to a file in tests/data/
fn test_data_path(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(filename)
}

/// Build a minimal valid GenBank record string.
///
/// The ORIGIN section is formatted with 9-digit right-justified positions,
/// space-separated 10-character groups of lowercase sequence, 60 characters
/// per line. The LOCUS line includes the correct bp count and topology.
fn genbank_record(name: &str, topology: &str, seq: &str) -> String {
    let bp = seq.len();
    let seq_lower = seq.to_lowercase();
    let mut result = format!(
        "LOCUS       {:<16} {} bp    DNA     {}   UNK \n",
        name, bp, topology
    );
    result.push_str(&format!("DEFINITION  {name} test sequence.\n"));
    result.push_str("ORIGIN\n");

    // Format ORIGIN lines: 9-digit position, then 6 groups of 10 chars (60 per line)
    let chars: Vec<char> = seq_lower.chars().collect();
    let mut pos = 0;
    while pos < chars.len() {
        result.push_str(&format!("{:>9}", pos + 1));
        // Up to 6 groups of 10 per line (60 chars total)
        for group in 0..6 {
            let start = pos + group * 10;
            if start >= chars.len() {
                break;
            }
            let end = (start + 10).min(chars.len());
            result.push(' ');
            for ch in &chars[start..end] {
                result.push(*ch);
            }
        }
        result.push('\n');
        pos += 60;
    }

    result.push_str("//\n");
    result
}

/// Write content to a temp file with the given suffix. Returns the `NamedTempFile`
/// handle (file stays alive while the handle exists).
fn write_temp_file(content: &str, suffix: &str) -> NamedTempFile {
    let mut file = NamedTempFile::with_suffix(suffix).expect("failed to create temp file");
    file.write_all(content.as_bytes())
        .expect("failed to write temp file");
    file.flush().expect("failed to flush temp file");
    file
}

/// Generate a deterministic DNA sequence of the given length.
/// Uses a repeating pattern `[A, C, G, T, A, T, G, C]` -- NOT random.
fn random_dna(len: usize) -> Vec<u8> {
    const PATTERN: &[u8] = b"ACGTATGC";
    (0..len).map(|i| PATTERN[i % PATTERN.len()]).collect()
}

/// Build a synthetic ~1800bp sequence with 16S primer binding sites.
///
/// Layout:
/// - 100bp filler
/// - Forward primer 27F binding site at position 100 (20bp, M resolved to A)
/// - ~1480bp filler
/// - RC of reverse primer 1492R at position 1600 (22bp, Y resolved to T)
/// - 100bp filler
///
/// Returns (record, expected_start=100, expected_end=1622) using 0-based half-open.
fn build_16s_fixture() -> (SequenceRecord, usize, usize) {
    // Forward primer 27F with M=A: AGAGTTTGATCATGGCTCAG (20bp)
    let fwd_binding: &[u8] = b"AGAGTTTGATCATGGCTCAG";
    // Reverse primer 1492R with Y=T: TACGGTTACCTTGTTACGACTT (22bp)
    // RC of that reverse primer goes on the + strand:
    let rev_primer_resolved: &[u8] = b"TACGGTTACCTTGTTACGACTT";
    let rc_rev = reverse_complement(rev_primer_resolved);
    // rc_rev = AAGTCGTAACAAGGTAACCGTA (22bp)

    let flank_left = random_dna(100);
    let filler_mid = random_dna(1480);
    let flank_right = random_dna(100);

    let mut seq = Vec::with_capacity(1722);
    seq.extend_from_slice(&flank_left); // 0..100
    seq.extend_from_slice(fwd_binding); // 100..120
    seq.extend_from_slice(&filler_mid); // 120..1600
    seq.extend_from_slice(&rc_rev); // 1600..1622
    seq.extend_from_slice(&flank_right); // 1622..1722

    let record = SequenceRecord {
        header: "synthetic_ecoli_16S".to_string(),
        sequence: seq,
        source_file: PathBuf::from("synthetic_16s.fasta"),
    };

    // Product spans from forward start (100) to RC-reverse end (1622)
    (record, 100, 1622)
}

/// Build a 500bp circular sequence for wrap-around product testing.
///
/// Layout:
/// - RC of reverse primer at position 30..44 (14bp)
/// - Forward primer at position 450..464 (14bp)
/// - All other positions filled with deterministic filler
///
/// Product wraps from 450 -> end(500) -> start -> 44 = ~94bp
///
/// Returns (record, forward_primer_bytes, reverse_primer_bytes).
fn build_circular_wrap_fixture() -> (SequenceRecord, Vec<u8>, Vec<u8>) {
    let fwd_primer = b"GCTAGCTAGCTAAC".to_vec(); // 14bp
    let rev_primer = b"AATTCCGGAATTCC".to_vec(); // 14bp
    let rc_rev = reverse_complement(&rev_primer);

    let mut seq = random_dna(500);

    // Place RC of reverse primer at position 30..44
    seq[30..44].copy_from_slice(&rc_rev);
    // Place forward primer at position 450..464
    seq[450..464].copy_from_slice(&fwd_primer);

    let record = SequenceRecord {
        header: "circular_plasmid_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("plasmid.gb"),
    };

    (record, fwd_primer, rev_primer)
}

/// Build a ~200bp sequence where primers only match on the RC strand.
///
/// Biology recap (per pcr.rs `find_rc_strand_products`):
/// - The entire sequence is reverse-complemented
/// - The same search logic applies: find forward primer, then RC(reverse) downstream
///
/// So on the RC strand, we need: forward primer at some position, then RC(reverse)
/// downstream. This means on the + strand we place:
/// - The reverse primer at positions 20..32 (on RC strand this becomes RC(reverse) near the end)
/// - RC(forward) at positions 140..152 (on RC strand this becomes the forward primer near the beginning)
///
/// On RC strand (200bp reversed):
/// - forward primer appears near position ~48 (200 - 152 = 48)
/// - RC(reverse) appears near position ~168 (200 - 32 = 168)
/// - Product spans ~48..~180 on the RC strand
///
/// Returns (record, forward_primer, reverse_primer, plus_strand_rc_fwd_start, plus_strand_rev_start).
fn build_rc_strand_fixture() -> (SequenceRecord, Vec<u8>, Vec<u8>, usize, usize) {
    let fwd_primer = b"CCTTAAGGCCTT".to_vec(); // 12bp
    let rev_primer = b"GGAATTCCGGAA".to_vec(); // 12bp
    let rc_fwd = reverse_complement(&fwd_primer);
    // rc_fwd = AAGGCCTTAAGG (12bp)

    let mut seq = random_dna(200);

    // Place reverse primer at + strand positions 20..32
    // On RC strand this becomes RC(reverse) near position 200-32 = 168
    seq[20..32].copy_from_slice(&rev_primer);

    // Place RC(forward) at + strand positions 140..152
    // On RC strand this becomes forward primer near position 200-152 = 48
    seq[140..152].copy_from_slice(&rc_fwd);

    let record = SequenceRecord {
        header: "rc_strand_test".to_string(),
        sequence: seq,
        source_file: PathBuf::from("rc_test.fasta"),
    };

    (record, fwd_primer, rev_primer, 140, 20)
}

// ============================================================================
// Placeholder test -- replaced in Task 2
// ============================================================================

#[test]
fn helpers_compile() {}
