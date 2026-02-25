use amplirust::matcher::PrimerMatch;
use amplirust::output::{write_fasta_record, write_tsv};
use amplirust::pcr::PcrProduct;
use insta::assert_snapshot;
use sassy::Strand;
use tempfile::NamedTempFile;

// ---------------------------------------------------------------------------
// Deterministic product construction helpers
// ---------------------------------------------------------------------------

/// Forward strand product with a short (20 bp) sequence.
fn make_fwd_product() -> PcrProduct {
    PcrProduct {
        reference_header: "NZ_CP012345.1 Escherichia coli strain K-12".to_string(),
        source_file: "test_genome.fasta".to_string(),
        primer_name: "27F_1492R".to_string(),
        sequence: b"AGAGTTTGATCATGGCTCAG".to_vec(),
        full_length: 20,
        fwd_match: PrimerMatch {
            start: 100,
            end: 120,
            edit_distance: 0,
            strand: Strand::Fwd,
            cigar: "20=".to_string(),
            identity: 1.0,
        },
        rev_match: PrimerMatch {
            start: 110,
            end: 120,
            edit_distance: 1,
            strand: Strand::Fwd,
            cigar: "9=1X".to_string(),
            identity: 0.9,
        },
        strand: Strand::Fwd,
        is_circular_wrap: false,
        original_start: 100,
        original_end: 120,
        case_number: 1,
    }
}

/// RC strand product (same data but with Rc strand markers).
fn make_rc_product() -> PcrProduct {
    PcrProduct {
        reference_header: "NZ_CP012345.1 Escherichia coli strain K-12".to_string(),
        source_file: "test_genome.fasta".to_string(),
        primer_name: "27F_1492R".to_string(),
        sequence: b"AGAGTTTGATCATGGCTCAG".to_vec(),
        full_length: 20,
        fwd_match: PrimerMatch {
            start: 100,
            end: 120,
            edit_distance: 0,
            strand: Strand::Rc,
            cigar: "20=".to_string(),
            identity: 1.0,
        },
        rev_match: PrimerMatch {
            start: 110,
            end: 120,
            edit_distance: 1,
            strand: Strand::Rc,
            cigar: "9=1X".to_string(),
            identity: 0.9,
        },
        strand: Strand::Rc,
        is_circular_wrap: false,
        original_start: 100,
        original_end: 120,
        case_number: 2,
    }
}

/// Circular wrap product (forward strand, wrapping around genome origin).
fn make_circular_wrap_product() -> PcrProduct {
    PcrProduct {
        reference_header: "NZ_CP012345.1 Escherichia coli strain K-12".to_string(),
        source_file: "test_genome.fasta".to_string(),
        primer_name: "27F_1492R".to_string(),
        sequence: b"AGAGTTTGATCATGGCTCAG".to_vec(),
        full_length: 20,
        fwd_match: PrimerMatch {
            start: 4990,
            end: 5010,
            edit_distance: 0,
            strand: Strand::Fwd,
            cigar: "20=".to_string(),
            identity: 1.0,
        },
        rev_match: PrimerMatch {
            start: 5000,
            end: 5010,
            edit_distance: 1,
            strand: Strand::Fwd,
            cigar: "9=1X".to_string(),
            identity: 0.9,
        },
        strand: Strand::Fwd,
        is_circular_wrap: true,
        original_start: 4990,
        original_end: 15,
        case_number: 1,
    }
}

/// Product with a 100 bp sequence to test FASTA line wrapping at 80 chars.
fn make_long_sequence_product() -> PcrProduct {
    PcrProduct {
        reference_header: "NZ_CP012345.1 Escherichia coli strain K-12".to_string(),
        source_file: "test_genome.fasta".to_string(),
        primer_name: "27F_1492R".to_string(),
        sequence: b"ACGT".repeat(25), // 100 bytes
        full_length: 100,
        fwd_match: PrimerMatch {
            start: 100,
            end: 200,
            edit_distance: 0,
            strand: Strand::Fwd,
            cigar: "100=".to_string(),
            identity: 1.0,
        },
        rev_match: PrimerMatch {
            start: 190,
            end: 200,
            edit_distance: 1,
            strand: Strand::Fwd,
            cigar: "9=1X".to_string(),
            identity: 0.9,
        },
        strand: Strand::Fwd,
        is_circular_wrap: false,
        original_start: 100,
        original_end: 200,
        case_number: 1,
    }
}

// ---------------------------------------------------------------------------
// FASTA snapshot tests (TEST-10)
// ---------------------------------------------------------------------------

#[test]
fn fasta_forward_strand() {
    let mut buf = Vec::new();
    write_fasta_record(&mut buf, &make_fwd_product()).unwrap();
    let content = String::from_utf8(buf).unwrap();
    assert_snapshot!("fasta_forward_strand", content);
}

#[test]
fn fasta_rc_strand() {
    let mut buf = Vec::new();
    write_fasta_record(&mut buf, &make_rc_product()).unwrap();
    let content = String::from_utf8(buf).unwrap();
    assert_snapshot!("fasta_rc_strand", content);
}

#[test]
fn fasta_circular_wrap() {
    let mut buf = Vec::new();
    write_fasta_record(&mut buf, &make_circular_wrap_product()).unwrap();
    let content = String::from_utf8(buf).unwrap();
    assert_snapshot!("fasta_circular_wrap", content);
}

#[test]
fn fasta_long_sequence_wrapping() {
    let mut buf = Vec::new();
    write_fasta_record(&mut buf, &make_long_sequence_product()).unwrap();
    let content = String::from_utf8(buf).unwrap();

    // Verify line wrapping at 80 chars
    let lines: Vec<&str> = content.lines().collect();
    assert!(
        lines.len() >= 3,
        "Should have header + at least 2 sequence lines"
    );
    for line in &lines[1..] {
        assert!(
            line.len() <= 80,
            "Sequence line should be max 80 chars, got {}",
            line.len()
        );
    }

    assert_snapshot!("fasta_long_sequence", content);
}

// ---------------------------------------------------------------------------
// TSV snapshot tests (TEST-11)
// ---------------------------------------------------------------------------

#[test]
fn tsv_output_with_header() {
    let products = vec![make_fwd_product()];
    let tmp = NamedTempFile::with_suffix(".tsv").unwrap();
    write_tsv(&products, tmp.path()).unwrap();
    let content = std::fs::read_to_string(tmp.path()).unwrap();

    // Sanity: header row should have 19 columns
    let header_cols: Vec<&str> = content.lines().next().unwrap().split('\t').collect();
    assert_eq!(header_cols.len(), 19, "TSV header should have 19 columns");

    assert_snapshot!("tsv_output_with_header", content);
}

#[test]
fn tsv_multi_product() {
    let products = vec![
        make_fwd_product(),
        make_rc_product(),
        make_circular_wrap_product(),
    ];
    let tmp = NamedTempFile::with_suffix(".tsv").unwrap();
    write_tsv(&products, tmp.path()).unwrap();
    let content = std::fs::read_to_string(tmp.path()).unwrap();

    // Should have header + 3 data rows
    let line_count = content.lines().count();
    assert_eq!(
        line_count, 4,
        "Should have 1 header + 3 data rows, got {line_count}"
    );

    // Each data row should have 19 columns
    for (i, line) in content.lines().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cols.len(),
            19,
            "Row {i} should have 19 columns, got {}",
            cols.len()
        );
    }

    assert_snapshot!("tsv_multi_product", content);
}
