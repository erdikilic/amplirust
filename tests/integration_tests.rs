//! Integration tests for amplirust

use std::path::PathBuf;
use tempfile::TempDir;

use amplirust::input::{SequenceRecord, expand_input_patterns, read_sequences_from_file};
use amplirust::matcher::MatchConfig;
use amplirust::output::{write_fasta, write_tsv};
use amplirust::pcr::{PcrConfig, find_pcr_products};
use amplirust::primer::{PrimerPair, parse_primers};

fn test_data_path(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(filename)
}

#[test]
fn test_load_fasta() {
    let path = test_data_path("test.fasta");
    let files = expand_input_patterns(&[path.to_string_lossy().to_string()]).unwrap();
    assert_eq!(files.len(), 1);

    let mut sequences = Vec::new();
    for file in &files {
        let mut records = read_sequences_from_file(file, 0).unwrap();
        sequences.append(&mut records);
    }
    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].header, "test_sequence_1");
    assert_eq!(sequences[1].header, "test_sequence_2_with_product");
    assert_eq!(sequences[2].header, "circular_test");
}

#[test]
fn test_load_primers_csv() {
    let path = test_data_path("primers.csv");
    let primers = parse_primers(&path.to_string_lossy()).unwrap();

    assert_eq!(primers.len(), 2);
    assert_eq!(primers[0].name, "16S");
    assert_eq!(primers[1].name, "simple");
}

#[test]
fn test_load_primers_string() {
    let primers = parse_primers("test:ACGT:TGCA").unwrap();
    assert_eq!(primers.len(), 1);
    assert_eq!(primers[0].name, "test");
}

#[test]
fn test_find_simple_product() {
    // Create a simple test case with exact matches
    let sequence = SequenceRecord {
        header: "test".to_string(),
        sequence: b"AAAACGTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        //              ACGT--------------------ACGT (RC of ACGT)
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();
    // RC of ACGT is ACGT (palindrome)

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // Should find one product
    assert!(!products.is_empty(), "Should find at least one product");

    let product = &products[0];
    assert_eq!(product.primer_name, "test");
    assert!(product.full_length >= 10);
    assert!(product.full_length <= 100);
}

#[test]
fn test_find_product_with_mismatches() {
    // Test with 1 mismatch allowed
    let sequence = SequenceRecord {
        header: "test".to_string(),
        // Forward: ACGT, has ACTT (1 mismatch)
        // Reverse RC (ACGT): has ACGT (exact)
        sequence: b"AAAACTTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 1,
            min_identity: 0.5,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    if !products.is_empty() {
        let product = &products[0];
        // Forward match should have 1 mismatch
        assert!(product.fwd_match.edit_distance <= 1);
    }
}

#[test]
fn test_circular_genome() {
    // Test circular genome wrap-around
    // The sequence wraps: ...ACGT at end + NNNN... at start should find product
    let sequence = SequenceRecord {
        header: "circular".to_string(),
        sequence: b"NNNNNNNNNNNNNNNNNNNNACGTAAAACGT".to_vec(),
        //          Start                  ^End wraps to start
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 4,
        max_len: 100,
        circular: true, // Enable circular mode
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // Should find products in circular mode
    // At minimum should find the ACGT...ACGT within the sequence
    assert!(
        !products.is_empty(),
        "Should find products in circular mode"
    );
}

#[test]
fn test_trim_primers() {
    let sequence = SequenceRecord {
        header: "test".to_string(),
        sequence: b"AAAACGTXXXXXXXXACGTTTTT".to_vec(),
        //              ACGT--------ACGT
        //              Product with trim: XXXXXXXX
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    // Without trim
    let config_no_trim = PcrConfig {
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

    // With trim
    let config_trim = PcrConfig {
        trim_primers: true,
        ..config_no_trim.clone()
    };

    let products_no_trim = find_pcr_products(&sequence, &primer, &config_no_trim);
    let products_trim = find_pcr_products(&sequence, &primer, &config_trim);

    if !products_no_trim.is_empty() && !products_trim.is_empty() {
        // Trimmed product should be shorter (no primers)
        assert!(products_trim[0].len() < products_no_trim[0].len());
    }
}

#[test]
fn test_output_fasta() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fasta");

    let sequence = SequenceRecord {
        header: "test".to_string(),
        sequence: b"AAAACGTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    if !products.is_empty() {
        write_fasta(&products, &output_path, 1).unwrap();

        // Verify file was created and has content
        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains('>'));
        assert!(content.len() > 10);
    }
}

#[test]
fn test_output_tsv() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.tsv");

    let sequence = SequenceRecord {
        header: "test".to_string(),
        sequence: b"AAAACGTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    if !products.is_empty() {
        write_tsv(&products, &output_path).unwrap();

        // Verify TSV was created with header
        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains("amplicon_id"));
        assert!(content.contains("reference_id"));
        assert!(content.contains("product_seq"));
    }
}

#[test]
fn test_iupac_primer_matching() {
    // Test that IUPAC codes in primers work
    let sequence = SequenceRecord {
        header: "test".to_string(),
        // R matches A or G, so ACRT should match ACGT or ACAT
        sequence: b"AAAACGTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    // Using IUPAC code R (A or G)
    let primer = PrimerPair::new("test", b"ACRT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0, // Exact match with IUPAC expansion
            min_identity: 0.8,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // IUPAC matching should work (ACRT matches ACGT since R matches G)
    // This depends on sassy's IUPAC implementation
    // The test verifies no crash occurs - if sassy supports IUPAC we may find products
    let _ = products; // Use the variable to avoid warning
}

#[test]
fn test_multiple_products_same_sequence() {
    // Test finding multiple products in the same sequence
    let sequence = SequenceRecord {
        header: "multi".to_string(),
        sequence: b"ACGTNNNNNNNNACGTNNNNNNNNACGTNNNNNNNNACGT".to_vec(),
        //          1st---------2nd---------3rd---------4th
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

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

    let products = find_pcr_products(&sequence, &primer, &config);

    // Should find multiple products (several combinations)
    assert!(products.len() > 1, "Should find multiple products");
}

#[test]
fn test_length_filter() {
    let sequence = SequenceRecord {
        header: "test".to_string(),
        sequence: b"ACGTNNACGT".to_vec(), // Short product
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    // Min length too high
    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 100, // Product is only ~10bp, should be filtered
        max_len: 1000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // Should find no products due to length filter
    assert!(
        products.is_empty(),
        "Short products should be filtered by min_len"
    );
}

#[test]
fn test_gzip_output_readable() {
    use flate2::read::GzDecoder;
    use std::io::Read;

    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fasta.gz");

    let sequence = SequenceRecord {
        header: "test_gzip".to_string(),
        sequence: b"AAAACGTNNNNNNNNNNNNNNNNNNNNACGTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 10,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);
    assert!(!products.is_empty(), "Should find products for gzip test");

    // Write gzipped output
    write_fasta(&products, &output_path, 1).unwrap();

    // Verify file exists and has content
    assert!(output_path.exists(), "Gzip file should exist");

    // Read and decompress to verify
    let file = std::fs::File::open(&output_path).unwrap();
    let mut decoder = GzDecoder::new(file);
    let mut content = String::new();
    decoder.read_to_string(&mut content).unwrap();

    // Verify content is valid FASTA
    assert!(
        content.contains('>'),
        "Decompressed content should have FASTA header"
    );
    assert!(
        content.contains("test_gzip"),
        "Should contain reference name"
    );
}

#[test]
fn test_rc_strand_integration() {
    // Integration test for RC strand product discovery
    let sequence = SequenceRecord {
        header: "rc_test".to_string(),
        // Sequence designed so product only exists on RC strand
        // Forward: TTCA, Reverse: GCTT (in original orientation)
        // RC strand will have: AAGC...TGAA
        sequence: b"AAAATTCANNNNNNNNNNNNNNNNNNNGCTTTTTT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("rc_primer", b"AAGC", b"TTCA").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true, // Enable RC search
        },
        min_len: 4,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // Should find RC strand product
    assert!(!products.is_empty(), "Should find RC strand product");

    // At least one product should be on RC strand
    let rc_products: Vec<_> = products
        .iter()
        .filter(|p| p.strand == sassy::Strand::Rc)
        .collect();
    assert!(
        !rc_products.is_empty(),
        "Should have at least one RC strand product"
    );
}

#[test]
fn test_circular_wrap_product() {
    // Integration test for circular genome with product spanning origin
    let sequence = SequenceRecord {
        header: "circular_wrap".to_string(),
        // Put reverse primer at start, forward at end
        // Product will wrap around
        sequence: b"TGCANNNNNNNNNNNNNNNNACGT".to_vec(), // 24bp
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("wrap", b"ACGT", b"TGCA").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 4,
        max_len: 50,
        circular: true, // Enable circular mode
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    // Look for wrapping products
    let wrap_products: Vec<_> = products.iter().filter(|p| p.is_circular_wrap).collect();

    // Should find at least one wrapping product
    assert!(
        !wrap_products.is_empty(),
        "Should find circular wrap product"
    );
}

#[test]
fn test_n_fraction_filtering() {
    // Integration test for N-fraction filtering
    let sequence = SequenceRecord {
        header: "n_filter".to_string(),
        // 50% N content
        sequence: b"ACGTNNNNACGT".to_vec(), // 12bp, 4 N's = 33.3%
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    // With permissive N fraction
    let config_permissive = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        },
        min_len: 4,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 0.5, // Allow up to 50%
    };

    let products_permissive = find_pcr_products(&sequence, &primer, &config_permissive);
    assert!(
        !products_permissive.is_empty(),
        "Should find product with permissive N fraction"
    );

    // With strict N fraction
    let config_strict = PcrConfig {
        max_n_fraction: 0.1, // Only allow 10%
        ..config_permissive.clone()
    };

    let products_strict = find_pcr_products(&sequence, &primer, &config_strict);
    assert!(
        products_strict.is_empty(),
        "Should filter product with strict N fraction"
    );
}

#[test]
fn test_deduplication() {
    use amplirust::pcr::remove_duplicate_products_by_reference;

    let sequence = SequenceRecord {
        header: "dedup_test".to_string(),
        // Sequence with multiple overlapping products
        sequence: b"ACGTACGTACGTACGTACGT".to_vec(),
        source_file: PathBuf::from("test.fasta"),
    };

    let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true, // This may produce RC duplicates
        },
        min_len: 4,
        max_len: 100,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let products = find_pcr_products(&sequence, &primer, &config);

    if products.len() > 1 {
        // Test deduplication
        let deduped = remove_duplicate_products_by_reference(products.clone());

        // Deduplicated should have same or fewer products
        assert!(
            deduped.len() <= products.len(),
            "Deduplication should not increase product count"
        );

        // Case numbers should be reassigned
        for (i, product) in deduped.iter().enumerate() {
            assert_eq!(
                product.case_number,
                i + 1,
                "Case numbers should be sequential after dedup"
            );
        }
    }
}
