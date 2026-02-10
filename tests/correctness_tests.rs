//! Biological correctness integration tests for amplirust.
//!
//! These tests prove the core biological correctness claim: every PCR product
//! reported is biologically correct with accurate coordinates and sequence content.
//!
//! TEST-01: 16S primer amplicon at expected coordinates
//! TEST-02: Circular plasmid wrap-around detection
//! TEST-03: Multi-record `GenBank` independent per-record results
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

/// Build a minimal valid `GenBank` record string.
///
/// The ORIGIN section is formatted with 9-digit right-justified positions,
/// space-separated 10-character groups of lowercase sequence, 60 characters
/// per line. The LOCUS line includes the correct bp count and topology.
fn genbank_record(name: &str, topology: &str, seq: &str) -> String {
    use std::fmt::Write as FmtWrite;

    let bp = seq.len();
    let seq_lower = seq.to_lowercase();
    let mut result = format!(
        "LOCUS       {name:<16} {bp} bp    DNA     {topology}   UNK \n"
    );
    let _ = writeln!(result, "DEFINITION  {name} test sequence.");
    result.push_str("ORIGIN\n");

    // Format ORIGIN lines: 9-digit position, then 6 groups of 10 chars (60 per line)
    let chars: Vec<char> = seq_lower.chars().collect();
    let mut pos = 0;
    while pos < chars.len() {
        let _ = write!(result, "{:>9}", pos + 1);
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
/// Returns `(record, expected_start=100, expected_end=1622)` using 0-based half-open.
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
/// Returns `(record, forward_primer_bytes, reverse_primer_bytes)`.
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
/// Returns `(record, forward_primer, reverse_primer, plus_strand_rc_fwd_start, plus_strand_rev_start)`.
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
// TEST-01: 16S primer amplicon at expected coordinates
// ============================================================================

#[test]
fn test_16s_primers_find_amplicon_at_expected_coordinates() {
    let (record, expected_start, expected_end) = build_16s_fixture();

    let primer = PrimerPair::new(
        "27F_1492R",
        b"AGAGTTTGATCMTGGCTCAG",
        b"TACGGYTACCTTGTTACGACTT",
    )
    .unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 1,
            min_identity: 0.8,
            search_rc: false,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    // UNCONDITIONAL: must find exactly one amplicon
    assert_eq!(
        products.len(),
        1,
        "Expected exactly 1 amplicon from 16S primer pair, got {}",
        products.len()
    );

    let product = &products[0];
    assert_eq!(
        product.original_start, expected_start,
        "Amplicon start should be {expected_start}, got {}",
        product.original_start
    );
    assert_eq!(
        product.original_end, expected_end,
        "Amplicon end should be {expected_end}, got {}",
        product.original_end
    );

    let expected_len = expected_end - expected_start;
    assert_eq!(
        product.full_length, expected_len,
        "Amplicon length should be {expected_len}, got {}",
        product.full_length
    );
    assert_eq!(
        product.strand,
        Strand::Fwd,
        "16S amplicon should be on forward strand"
    );
}

#[test]
fn test_16s_primers_from_csv_match_same_as_inline() {
    let (record, expected_start, expected_end) = build_16s_fixture();

    // Load primers from CSV file
    let primers = parse_primers(&test_data_path("primers_16s.csv").to_string_lossy()).unwrap();
    assert_eq!(primers.len(), 1, "CSV should contain exactly 1 primer pair");

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 1,
            min_identity: 0.8,
            search_rc: false,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primers[0], &config);

    // UNCONDITIONAL: CSV primers should produce the same result as inline
    assert_eq!(
        products.len(),
        1,
        "CSV primers should find exactly 1 amplicon, got {}",
        products.len()
    );
    assert_eq!(products[0].original_start, expected_start);
    assert_eq!(products[0].original_end, expected_end);
}

#[test]
fn test_16s_amplicon_sequence_starts_with_forward_primer() {
    let (record, _expected_start, _expected_end) = build_16s_fixture();

    let primer = PrimerPair::new(
        "27F_1492R",
        b"AGAGTTTGATCMTGGCTCAG",
        b"TACGGYTACCTTGTTACGACTT",
    )
    .unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 1,
            min_identity: 0.8,
            search_rc: false,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);
    assert_eq!(products.len(), 1, "Expected exactly 1 amplicon");

    let product = &products[0];

    // Forward primer binding site (M resolved to A): AGAGTTTGATCATGGCTCAG
    let fwd_binding = b"AGAGTTTGATCATGGCTCAG";
    assert!(
        product.sequence.starts_with(fwd_binding),
        "Product sequence should start with forward primer binding site.\n\
         Expected start: {:?}\n\
         Actual start:   {:?}",
        String::from_utf8_lossy(fwd_binding),
        String::from_utf8_lossy(&product.sequence[..20.min(product.sequence.len())]),
    );

    // RC of reverse primer (Y=T): AAGTCGTAACAAGGTAACCGTA (22bp)
    let rev_primer_resolved = b"TACGGTTACCTTGTTACGACTT";
    let rc_rev = reverse_complement(rev_primer_resolved);
    assert!(
        product.sequence.ends_with(&rc_rev),
        "Product sequence should end with RC of reverse primer binding site.\n\
         Expected end: {:?}\n\
         Actual end:   {:?}",
        String::from_utf8_lossy(&rc_rev),
        String::from_utf8_lossy(
            &product.sequence[product.sequence.len().saturating_sub(22)..]
        ),
    );
}

// ============================================================================
// TEST-02: Circular plasmid wrap-around
// ============================================================================

#[test]
fn test_circular_wrap_product_found() {
    let (record, fwd_primer, rev_primer) = build_circular_wrap_fixture();

    let primer = PrimerPair::new("circ_test", &fwd_primer, &rev_primer).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 200,
        circular: true,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    // UNCONDITIONAL: at least one product must have is_circular_wrap = true
    let wrap_products: Vec<_> = products.iter().filter(|p| p.is_circular_wrap).collect();
    assert!(
        !wrap_products.is_empty(),
        "Expected at least one circular wrap product. Found {} total products, none wrapping.\n\
         Products: {:?}",
        products.len(),
        products
            .iter()
            .map(|p| (p.original_start, p.original_end, p.is_circular_wrap))
            .collect::<Vec<_>>()
    );
}

#[test]
fn test_circular_wrap_product_sequence_content() {
    let (record, fwd_primer, rev_primer) = build_circular_wrap_fixture();
    let rc_rev = reverse_complement(&rev_primer);

    let primer = PrimerPair::new("circ_test", &fwd_primer, &rev_primer).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 200,
        circular: true,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);
    let wrap_products: Vec<_> = products.iter().filter(|p| p.is_circular_wrap).collect();
    assert!(
        !wrap_products.is_empty(),
        "Expected at least one circular wrap product"
    );

    let product = &wrap_products[0];

    // The wrapping product sequence should start with the forward primer
    assert!(
        product.sequence.starts_with(&fwd_primer),
        "Wrap product should start with forward primer.\n\
         Expected: {:?}\n\
         Got:      {:?}",
        String::from_utf8_lossy(&fwd_primer),
        String::from_utf8_lossy(
            &product.sequence[..fwd_primer.len().min(product.sequence.len())]
        ),
    );

    // The wrapping product sequence should end with RC(reverse)
    assert!(
        product.sequence.ends_with(&rc_rev),
        "Wrap product should end with RC of reverse primer.\n\
         Expected: {:?}\n\
         Got:      {:?}",
        String::from_utf8_lossy(&rc_rev),
        String::from_utf8_lossy(
            &product.sequence[product.sequence.len().saturating_sub(rc_rev.len())..]
        ),
    );

    // Product length should be: from fwd (450) to end (500) + start to rc_rev end (44)
    // = 50 + 44 = 94bp
    let expected_len = (500 - 450) + 44;
    assert_eq!(
        product.sequence.len(),
        expected_len,
        "Wrap product length should be {expected_len}, got {}",
        product.sequence.len()
    );
}

// ============================================================================
// TEST-03: Multi-record GenBank independent per-record results
// ============================================================================

#[test]
fn test_multi_record_genbank_produces_per_record_results() {
    // Use simple 14bp ACGT-only primers for this test
    let fwd_primer = b"GCTAGCTAGCTAAC";
    let rev_primer = b"AATTCCGGAATTCC";
    let rc_rev = reverse_complement(rev_primer);

    // Build Record 1: linear, 300bp, has primer sites
    let mut seq1_bytes = random_dna(300);
    seq1_bytes[50..64].copy_from_slice(fwd_primer); // Forward at 50
    seq1_bytes[200..214].copy_from_slice(&rc_rev); // RC(reverse) at 200
    let seq1: String = seq1_bytes.iter().map(|&b| b as char).collect();

    // Build Record 2: circular, 400bp, has primer sites at different positions
    let mut seq2_bytes = random_dna(400);
    seq2_bytes[100..114].copy_from_slice(fwd_primer); // Forward at 100
    seq2_bytes[300..314].copy_from_slice(&rc_rev); // RC(reverse) at 300
    let seq2: String = seq2_bytes.iter().map(|&b| b as char).collect();

    // Build Record 3: linear, 200bp, no primer sites (all filler)
    let seq3_bytes = random_dna(200);
    let seq3: String = seq3_bytes.iter().map(|&b| b as char).collect();

    // Build the multi-record GenBank file
    let mut gb_content = String::new();
    gb_content.push_str(&genbank_record("Rec1_linear", "linear", &seq1));
    gb_content.push_str(&genbank_record("Rec2_circular", "circular", &seq2));
    gb_content.push_str(&genbank_record("Rec3_no_match", "linear", &seq3));

    let temp_file = write_temp_file(&gb_content, ".gb");

    // Read records from the GenBank file
    let records = read_sequences_from_file(temp_file.path(), 0).unwrap();

    // UNCONDITIONAL: should have exactly 3 records
    assert_eq!(
        records.len(),
        3,
        "Expected 3 records from multi-record GenBank, got {}",
        records.len()
    );

    let primer_pair = PrimerPair::new("multi_test", fwd_primer, rev_primer).unwrap();

    let config_linear = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let config_circular = PcrConfig {
        circular: true,
        ..config_linear.clone()
    };

    // Record 1: linear, should find exactly 1 product
    let products1 = find_pcr_products(&records[0], &primer_pair, &config_linear);
    assert_eq!(
        products1.len(),
        1,
        "Record 1 (Rec1_linear) should produce exactly 1 product, got {}",
        products1.len()
    );

    // Record 2: circular, should find at least 1 product
    let products2 = find_pcr_products(&records[1], &primer_pair, &config_circular);
    assert!(
        !products2.is_empty(),
        "Record 2 (Rec2_circular) should produce at least 1 product, got 0"
    );

    // Record 3: no primers, should find 0 products
    let products3 = find_pcr_products(&records[2], &primer_pair, &config_linear);
    assert_eq!(
        products3.len(),
        0,
        "Record 3 (Rec3_no_match) should produce 0 products, got {}",
        products3.len()
    );

    // Verify Record 1 product has correct reference_header
    assert_eq!(
        products1[0].reference_header, "Rec1_linear",
        "Record 1 product should reference 'Rec1_linear', got '{}'",
        products1[0].reference_header
    );
}

// ============================================================================
// TEST-04: RC strand products
// ============================================================================

#[test]
fn test_rc_strand_product_has_correct_strand() {
    let (record, fwd_primer, rev_primer, _, _) = build_rc_strand_fixture();

    let primer = PrimerPair::new("rc_test", &fwd_primer, &rev_primer).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true,
        },
        min_len: 10,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);

    // UNCONDITIONAL: must find at least one product on RC strand
    let rc_products: Vec<_> = products.iter().filter(|p| p.strand == Strand::Rc).collect();
    assert!(
        !rc_products.is_empty(),
        "Expected at least one RC-strand product. Found {} total products, none on RC.\n\
         Products: {:?}",
        products.len(),
        products
            .iter()
            .map(|p| (p.strand, p.original_start, p.original_end))
            .collect::<Vec<_>>()
    );

    // Verify the strand field is correct
    for product in &rc_products {
        assert_eq!(product.strand, Strand::Rc);
    }
}

#[test]
fn test_rc_strand_product_sequence_is_rc_of_original_region() {
    let (record, fwd_primer, rev_primer, _, _) = build_rc_strand_fixture();

    let primer = PrimerPair::new("rc_test", &fwd_primer, &rev_primer).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true,
        },
        min_len: 10,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);
    let rc_products: Vec<_> = products.iter().filter(|p| p.strand == Strand::Rc).collect();
    assert!(
        !rc_products.is_empty(),
        "Expected at least one RC-strand product"
    );

    let product = &rc_products[0];

    // Extract the original + strand region
    let original_region = &record.sequence[product.original_start..product.original_end];
    let expected_sequence = reverse_complement(original_region);

    assert_eq!(
        product.sequence, expected_sequence,
        "RC product sequence should be RC of the original + strand region.\n\
         original_start={}, original_end={}\n\
         Product len: {}, Expected len: {}",
        product.original_start,
        product.original_end,
        product.sequence.len(),
        expected_sequence.len()
    );
}

#[test]
fn test_rc_strand_coordinates_map_to_plus_strand() {
    let (record, fwd_primer, rev_primer, _, _) = build_rc_strand_fixture();

    let primer = PrimerPair::new("rc_test", &fwd_primer, &rev_primer).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true,
        },
        min_len: 10,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&record, &primer, &config);
    let rc_products: Vec<_> = products.iter().filter(|p| p.strand == Strand::Rc).collect();
    assert!(
        !rc_products.is_empty(),
        "Expected at least one RC-strand product"
    );

    for product in &rc_products {
        // Coordinates should be valid: start < end, end <= sequence length
        assert!(
            product.original_start < product.original_end,
            "original_start ({}) should be less than original_end ({})",
            product.original_start,
            product.original_end
        );
        assert!(
            product.original_end <= record.sequence.len(),
            "original_end ({}) should be <= sequence length ({})",
            product.original_end,
            record.sequence.len()
        );
    }
}
