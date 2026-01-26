use crate::input::SequenceRecord;
use crate::matcher::{MatchConfig, PrimerMatch, PrimerMatcher};
use crate::primer::PrimerPair;
use crate::utils::{
    circular_to_original_pos, is_circular_wrap, make_circular_searchable, reverse_complement,
};
use sassy::Strand;

/// Configuration for PCR product detection
#[derive(Debug, Clone)]
pub struct PcrConfig {
    /// Matching configuration
    pub match_config: MatchConfig,
    /// Minimum product length (including primers)
    pub min_len: usize,
    /// Maximum product length (including primers)
    pub max_len: usize,
    /// Whether to treat sequences as circular
    pub circular: bool,
    /// Whether to trim primers from output
    pub trim_primers: bool,
    /// Maximum fraction of N bases allowed in product sequence
    pub max_n_fraction: f64,
}

impl Default for PcrConfig {
    fn default() -> Self {
        Self {
            match_config: MatchConfig::default(),
            min_len: 50,
            max_len: 5000,
            circular: false,
            trim_primers: false,
            max_n_fraction: 1.0,
        }
    }
}

/// A PCR product found by matching primers
#[derive(Debug, Clone)]
pub struct PcrProduct {
    /// Reference sequence header
    pub reference_header: String,
    /// Source file path
    pub source_file: String,
    /// Primer pair name
    pub primer_name: String,
    /// Product sequence (may be trimmed if trim_primers is set)
    pub sequence: Vec<u8>,
    /// Full product length (including primers, before trimming)
    pub full_length: usize,
    /// Forward primer match details
    pub fwd_match: PrimerMatch,
    /// Reverse primer match details
    pub rev_match: PrimerMatch,
    /// Strand (Fwd = + strand, Rc = - strand)
    pub strand: Strand,
    /// Whether product wraps around (circular genome)
    pub is_circular_wrap: bool,
    /// Start position in original sequence (0-based)
    pub original_start: usize,
    /// End position in original sequence (exclusive)
    pub original_end: usize,
    /// Case number for this reference header
    pub case_number: usize,
}

impl PcrProduct {
    /// Reference identifier (first whitespace-delimited token of the header)
    pub fn reference_id(&self) -> &str {
        self.reference_header
            .split_whitespace()
            .next()
            .unwrap_or(self.reference_header.as_str())
    }

    /// Generate the output header for this product
    pub fn header(&self) -> String {
        let strand_suffix = match self.strand {
            Strand::Fwd => "",
            Strand::Rc => "_rc",
        };
        let wrap_suffix = if self.is_circular_wrap { "_wrap" } else { "" };

        format!(
            "{}:{}{}{}:{}",
            self.reference_id(),
            self.primer_name,
            strand_suffix,
            wrap_suffix,
            self.case_number
        )
    }

    /// Get the product length
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if product is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Find all PCR products for a single sequence and primer pair
pub fn find_pcr_products(
    record: &SequenceRecord,
    primer: &PrimerPair,
    config: &PcrConfig,
) -> Vec<PcrProduct> {
    let mut products = Vec::new();
    let original_len = record.sequence.len();

    if original_len == 0 {
        return products;
    }

    // Prepare sequence (extend for circular if needed)
    let search_seq = if config.circular {
        make_circular_searchable(&record.sequence, config.max_len)
    } else {
        record.sequence.clone()
    };

    // Create matcher for forward strand only
    let fwd_config = MatchConfig {
        search_rc: false,
        ..config.match_config.clone()
    };
    let mut matcher = PrimerMatcher::new(fwd_config);

    // Find forward primer matches
    let fwd_matches = matcher.find_matches(&primer.forward, &search_seq);

    // For each forward match, look for reverse primer downstream
    // The reverse primer binds to the opposite strand, so we need to search
    // for its reverse complement in the sequence
    let rev_primer_rc = reverse_complement(&primer.reverse);

    for fwd_match in &fwd_matches {
        // Calculate the search region for reverse primer
        // It must be downstream of the forward primer
        let search_start = fwd_match.end;
        let search_end = (fwd_match.start + config.max_len).min(search_seq.len());

        if search_start >= search_end {
            continue;
        }

        let downstream_seq = &search_seq[search_start..search_end];

        // Create a new matcher for downstream search (always forward only here)
        let downstream_config = MatchConfig {
            search_rc: false, // We're searching for RC of reverse primer directly
            ..config.match_config.clone()
        };
        let mut downstream_matcher = PrimerMatcher::new(downstream_config);

        let rev_matches = downstream_matcher.find_matches(&rev_primer_rc, downstream_seq);

        for rev_match in rev_matches {
            // Adjust reverse match positions to original sequence coordinates
            let rev_start_abs = search_start + rev_match.start;
            let rev_end_abs = search_start + rev_match.end;

            // Calculate product boundaries
            let product_start = fwd_match.start;
            let product_end = rev_end_abs;
            let product_len = product_end - product_start;

            // Check length constraints
            if product_len < config.min_len || product_len > config.max_len {
                continue;
            }

            // Check for circular wrap
            let is_wrap = config.circular && is_circular_wrap(product_end - 1, original_len);

            // Extract product sequence
            let (sequence, original_start, original_end) = if config.trim_primers {
                // Trimmed: region between primers
                let trim_start = fwd_match.end;
                let trim_end = rev_start_abs;
                if trim_start >= trim_end {
                    continue; // No sequence between primers
                }
                let seq = search_seq[trim_start..trim_end].to_vec();
                (
                    seq,
                    circular_to_original_pos(trim_start, original_len),
                    circular_to_original_pos(trim_end, original_len),
                )
            } else {
                // Full product including primers
                let seq = search_seq[product_start..product_end].to_vec();
                (
                    seq,
                    circular_to_original_pos(product_start, original_len),
                    circular_to_original_pos(product_end, original_len),
                )
            };

            if n_fraction(&sequence) > config.max_n_fraction {
                continue;
            }

            // Adjust reverse match for storage (with absolute positions)
            let adjusted_rev_match = PrimerMatch {
                start: rev_start_abs,
                end: rev_end_abs,
                edit_distance: rev_match.edit_distance,
                strand: rev_match.strand,
                cigar: rev_match.cigar.clone(),
                identity: rev_match.identity,
            };

            products.push(PcrProduct {
                reference_header: record.header.clone(),
                source_file: record.source_file.display().to_string(),
                primer_name: primer.name.clone(),
                sequence,
                full_length: product_len,
                fwd_match: fwd_match.clone(),
                rev_match: adjusted_rev_match,
                strand: fwd_match.strand,
                is_circular_wrap: is_wrap,
                original_start,
                original_end,
                case_number: 0, // Will be assigned later
            });
        }
    }

    // If searching RC mode, also look for products on the reverse strand
    // This means: forward primer on RC strand, reverse primer on forward strand
    if config.match_config.search_rc {
        let rc_products =
            find_rc_strand_products(record, primer, config, &search_seq, original_len);
        products.extend(rc_products);
    }

    // Assign case numbers
    for (i, product) in products.iter_mut().enumerate() {
        product.case_number = i + 1;
    }

    products
}

/// Find PCR products on the reverse complement strand
fn find_rc_strand_products(
    record: &SequenceRecord,
    primer: &PrimerPair,
    config: &PcrConfig,
    search_seq: &[u8],
    original_len: usize,
) -> Vec<PcrProduct> {
    let mut products = Vec::new();

    // For RC strand products:
    // - Forward primer matches on RC strand (search for RC of forward primer)
    // - Reverse primer matches on forward strand downstream

    // Create forward-only matcher
    let fwd_config = MatchConfig {
        search_rc: false,
        ..config.match_config.clone()
    };
    let mut matcher = PrimerMatcher::new(fwd_config);

    // For reverse primer on + strand, we search for the reverse primer directly
    // (not its RC, because on the + strand it would bind as-is on the - strand template)
    // Actually, let's think about this more carefully:
    //
    // In PCR on - strand:
    // - Forward primer binds to + strand template (so we search its RC on sequence)
    // - Reverse primer binds to - strand template (so we search it directly)
    //
    // Wait, let me reconsider the biology:
    // For a product on the - strand:
    // - The forward primer anneals to the + strand (3'->5') and extends on - strand
    // - The reverse primer anneals to the - strand (3'->5') and extends on + strand
    //
    // If our sequence is the + strand:
    // - Forward primer RC will match where forward primer binds on + strand
    // - Reverse primer will match where reverse primer binds on - strand...
    //   but we're searching the + strand, so we need RC of reverse primer
    //
    // This is getting complex. Let's simplify:
    // For RC strand products, we reverse complement the entire sequence and
    // run the same algorithm. The positions are then adjusted.

    let rc_sequence = reverse_complement(search_seq);
    let rc_search_seq = if config.circular {
        make_circular_searchable(&rc_sequence[..original_len], config.max_len)
    } else {
        rc_sequence
    };

    // Now search on RC sequence with same logic as forward
    let fwd_matches = matcher.find_matches(&primer.forward, &rc_search_seq);
    let rev_primer_rc = reverse_complement(&primer.reverse);

    for fwd_match in &fwd_matches {
        let search_start = fwd_match.end;
        let search_end = (fwd_match.start + config.max_len).min(rc_search_seq.len());

        if search_start >= search_end {
            continue;
        }

        let downstream_seq = &rc_search_seq[search_start..search_end];
        let rev_matches = matcher.find_matches(&rev_primer_rc, downstream_seq);

        for rev_match in rev_matches {
            let rev_start_abs = search_start + rev_match.start;
            let rev_end_abs = search_start + rev_match.end;

            let product_start = fwd_match.start;
            let product_end = rev_end_abs;
            let product_len = product_end - product_start;

            if product_len < config.min_len || product_len > config.max_len {
                continue;
            }

            let is_wrap = config.circular && is_circular_wrap(product_end - 1, original_len);

            let (sequence, original_start, original_end) = if config.trim_primers {
                let trim_start = fwd_match.end;
                let trim_end = rev_start_abs;
                if trim_start >= trim_end {
                    continue;
                }
                let seq = rc_search_seq[trim_start..trim_end].to_vec();
                // Convert RC positions back to original coordinates
                let orig_end =
                    original_len.saturating_sub(circular_to_original_pos(trim_start, original_len));
                let orig_start =
                    original_len.saturating_sub(circular_to_original_pos(trim_end, original_len));
                (seq, orig_start, orig_end)
            } else {
                let seq = rc_search_seq[product_start..product_end].to_vec();
                let orig_end = original_len
                    .saturating_sub(circular_to_original_pos(product_start, original_len));
                let orig_start = original_len
                    .saturating_sub(circular_to_original_pos(product_end, original_len));
                (seq, orig_start.min(orig_end), orig_start.max(orig_end))
            };

            if n_fraction(&sequence) > config.max_n_fraction {
                continue;
            }

            let adjusted_rev_match = PrimerMatch {
                start: rev_start_abs,
                end: rev_end_abs,
                edit_distance: rev_match.edit_distance,
                strand: Strand::Rc, // Mark as from RC strand
                cigar: rev_match.cigar.clone(),
                identity: rev_match.identity,
            };

            let adjusted_fwd_match = PrimerMatch {
                strand: Strand::Rc,
                ..fwd_match.clone()
            };

            products.push(PcrProduct {
                reference_header: record.header.clone(),
                source_file: record.source_file.display().to_string(),
                primer_name: primer.name.clone(),
                sequence,
                full_length: product_len,
                fwd_match: adjusted_fwd_match,
                rev_match: adjusted_rev_match,
                strand: Strand::Rc,
                is_circular_wrap: is_wrap,
                original_start,
                original_end,
                case_number: 0,
            });
        }
    }

    products
}

/// Remove duplicate product sequences per reference header
pub fn remove_duplicate_products_by_reference(products: Vec<PcrProduct>) -> Vec<PcrProduct> {
    use std::collections::{HashMap, HashSet};

    let mut seen: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();
    let mut deduped = Vec::with_capacity(products.len());

    for product in products {
        let entry = seen.entry(product.reference_header.clone()).or_default();
        let canonical = canonical_sequence(&product.sequence);
        if entry.contains(&canonical) {
            continue;
        }
        entry.insert(canonical);
        deduped.push(product);
    }

    assign_case_numbers_by_reference(&mut deduped);
    deduped
}

fn assign_case_numbers_by_reference(products: &mut [PcrProduct]) {
    use std::collections::HashMap;

    let mut ref_counts: HashMap<&str, usize> = HashMap::new();
    for product in products.iter_mut() {
        let counter = ref_counts
            .entry(product.reference_header.as_str())
            .or_insert(0);
        *counter += 1;
        product.case_number = *counter;
    }
}

pub fn canonical_sequence(sequence: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(sequence);
    if rc.as_slice() < sequence {
        rc
    } else {
        sequence.to_vec()
    }
}

fn n_fraction(sequence: &[u8]) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    let n_count = sequence.iter().filter(|&&b| b == b'N').count();
    n_count as f64 / sequence.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn make_record(header: &str, seq: &[u8]) -> SequenceRecord {
        SequenceRecord {
            header: header.to_string(),
            sequence: seq.to_vec(),
            source_file: PathBuf::from("test.fasta"),
        }
    }

    #[test]
    fn test_simple_pcr_product() {
        let record = make_record("test", b"AAAACGTACGTACGTACGTTTTTT");
        //                                   ACGT............ACGT (RC=ACGT)
        let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

        let config = PcrConfig {
            match_config: MatchConfig {
                max_errors: 0,
                min_identity: 1.0,
                search_rc: false,
            },
            min_len: 8,
            max_len: 100,
            circular: false,
            trim_primers: false,
            max_n_fraction: 1.0,
        };

        let products = find_pcr_products(&record, &primer, &config);
        // Should find product from first ACGT to last ACGT
        assert!(!products.is_empty());
    }

    #[test]
    fn test_product_header() {
        let product = PcrProduct {
            reference_header: "chr1".to_string(),
            source_file: "test.fa".to_string(),
            primer_name: "16S".to_string(),
            sequence: b"ACGT".to_vec(),
            full_length: 100,
            fwd_match: PrimerMatch {
                start: 0,
                end: 4,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            rev_match: PrimerMatch {
                start: 96,
                end: 100,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            strand: Strand::Fwd,
            is_circular_wrap: false,
            original_start: 0,
            original_end: 100,
            case_number: 1,
        };

        assert_eq!(product.header(), "chr1:16S:1");
    }

    #[test]
    fn test_product_header_uses_reference_id() {
        let product = PcrProduct {
            reference_header: "NZ_CP172019.1 Bifidobacterium adolescentis".to_string(),
            source_file: "test.fa".to_string(),
            primer_name: "16S".to_string(),
            sequence: b"ACGT".to_vec(),
            full_length: 100,
            fwd_match: PrimerMatch {
                start: 0,
                end: 4,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            rev_match: PrimerMatch {
                start: 96,
                end: 100,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            strand: Strand::Fwd,
            is_circular_wrap: false,
            original_start: 0,
            original_end: 100,
            case_number: 1,
        };

        assert_eq!(product.header(), "NZ_CP172019.1:16S:1");
    }

    #[test]
    fn test_rc_product_header() {
        let product = PcrProduct {
            reference_header: "chr1".to_string(),
            source_file: "test.fa".to_string(),
            primer_name: "16S".to_string(),
            sequence: b"ACGT".to_vec(),
            full_length: 100,
            fwd_match: PrimerMatch {
                start: 0,
                end: 4,
                edit_distance: 0,
                strand: Strand::Rc,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            rev_match: PrimerMatch {
                start: 96,
                end: 100,
                edit_distance: 0,
                strand: Strand::Rc,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            strand: Strand::Rc,
            is_circular_wrap: false,
            original_start: 0,
            original_end: 100,
            case_number: 2,
        };

        assert_eq!(product.header(), "chr1:16S_rc:2");
    }

    #[test]
    fn test_remove_duplicate_products_by_reference() {
        let base = PcrProduct {
            reference_header: "chr1".to_string(),
            source_file: "test.fa".to_string(),
            primer_name: "16S".to_string(),
            sequence: b"ACGTTGCA".to_vec(),
            full_length: 8,
            fwd_match: PrimerMatch {
                start: 0,
                end: 4,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            rev_match: PrimerMatch {
                start: 96,
                end: 100,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            strand: Strand::Fwd,
            is_circular_wrap: false,
            original_start: 0,
            original_end: 8,
            case_number: 1,
        };

        let mut dup = base.clone();
        dup.strand = Strand::Rc;
        dup.sequence = reverse_complement(&base.sequence);
        let mut other_ref = base.clone();
        other_ref.reference_header = "chr2".to_string();

        let products = vec![base, dup, other_ref];
        let deduped = remove_duplicate_products_by_reference(products);

        assert_eq!(deduped.len(), 2);
        let chr1 = deduped
            .iter()
            .find(|p| p.reference_header == "chr1")
            .unwrap();
        let chr2 = deduped
            .iter()
            .find(|p| p.reference_header == "chr2")
            .unwrap();
        assert_eq!(chr1.case_number, 1);
        assert_eq!(chr2.case_number, 1);
    }

    #[test]
    fn test_max_n_fraction_filter() {
        let record = make_record("test", b"AAAACGTNNNNNNNNACGTTTTT");
        let primer = PrimerPair::new("test", b"ACGT", b"ACGT").unwrap();

        let config = PcrConfig {
            match_config: MatchConfig {
                max_errors: 0,
                min_identity: 1.0,
                search_rc: false,
            },
            min_len: 8,
            max_len: 100,
            circular: false,
            trim_primers: false,
            max_n_fraction: 0.2,
        };

        let products = find_pcr_products(&record, &primer, &config);
        assert!(products.is_empty());
    }

    #[test]
    fn test_rc_products_not_reversed() {
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AAAA");
        seq.extend_from_slice(b"AAGC");
        seq.extend_from_slice(b"CCCC");
        seq.extend_from_slice(b"TGAA");
        seq.extend_from_slice(b"GGGG");
        seq.extend_from_slice(b"TTCA");
        seq.extend_from_slice(b"TTTT");
        seq.extend_from_slice(b"GCTT");
        seq.extend_from_slice(b"CCC");

        let record = make_record("test", &seq);
        let primer = PrimerPair::new("test", b"AAGC", b"TTCA").unwrap();

        let config = PcrConfig {
            match_config: MatchConfig {
                max_errors: 0,
                min_identity: 1.0,
                search_rc: true,
            },
            min_len: 4,
            max_len: 200,
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

        for product in rc_products {
            let slice = &record.sequence[product.original_start..product.original_end];
            let expected = reverse_complement(slice);
            assert_eq!(product.sequence, expected);
        }
    }
}
