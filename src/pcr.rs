use sassy::Strand;
use crate::input::SequenceRecord;
use crate::matcher::{MatchConfig, PrimerMatch, PrimerMatcher};
use crate::primer::PrimerPair;
use crate::utils::{make_circular_searchable, reverse_complement, is_circular_wrap, circular_to_original_pos};

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
}

impl Default for PcrConfig {
    fn default() -> Self {
        Self {
            match_config: MatchConfig::default(),
            min_len: 50,
            max_len: 5000,
            circular: false,
            trim_primers: false,
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
    /// Case number for this reference+primer combination
    pub case_number: usize,
}

impl PcrProduct {
    /// Generate the output header for this product
    pub fn header(&self) -> String {
        let strand_suffix = match self.strand {
            Strand::Fwd => "",
            Strand::Rc => "_rc",
        };
        let wrap_suffix = if self.is_circular_wrap { "_wrap" } else { "" };
        
        format!(
            "{}_amplicon:{}{}{}:{}",
            self.reference_header,
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

    // Create matcher
    let mut matcher = PrimerMatcher::new(config.match_config.clone());

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
            search_rc: false,  // We're searching for RC of reverse primer directly
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
                    continue;  // No sequence between primers
                }
                let seq = search_seq[trim_start..trim_end].to_vec();
                (seq, circular_to_original_pos(trim_start, original_len), 
                 circular_to_original_pos(trim_end, original_len))
            } else {
                // Full product including primers
                let seq = search_seq[product_start..product_end].to_vec();
                (seq, circular_to_original_pos(product_start, original_len),
                 circular_to_original_pos(product_end, original_len))
            };

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
                case_number: 0,  // Will be assigned later
            });
        }
    }

    // If searching RC mode, also look for products on the reverse strand
    // This means: forward primer on RC strand, reverse primer on forward strand
    if config.match_config.search_rc {
        let rc_products = find_rc_strand_products(record, primer, config, &search_seq, original_len);
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
                let orig_end = original_len.saturating_sub(circular_to_original_pos(trim_start, original_len));
                let orig_start = original_len.saturating_sub(circular_to_original_pos(trim_end, original_len));
                (seq, orig_start, orig_end)
            } else {
                let seq = rc_search_seq[product_start..product_end].to_vec();
                let orig_end = original_len.saturating_sub(circular_to_original_pos(product_start, original_len));
                let orig_start = original_len.saturating_sub(circular_to_original_pos(product_end, original_len));
                (seq, orig_start.min(orig_end), orig_start.max(orig_end))
            };

            let adjusted_rev_match = PrimerMatch {
                start: rev_start_abs,
                end: rev_end_abs,
                edit_distance: rev_match.edit_distance,
                strand: Strand::Rc,  // Mark as from RC strand
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

/// Find all PCR products for multiple sequences and primers
pub fn find_all_products(
    records: &[SequenceRecord],
    primers: &[PrimerPair],
    config: &PcrConfig,
    show_progress: bool,
) -> Vec<PcrProduct> {
    use indicatif::{ProgressBar, ProgressStyle};
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    // Create progress bar if requested
    let pb = if show_progress && records.len() > 1 {
        let pb = ProgressBar::new(records.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} sequences ({eta})")
                .unwrap()
                .progress_chars("#>-"),
        );
        pb.set_message("Searching");
        Some(pb)
    } else {
        None
    };

    let completed = AtomicUsize::new(0);

    // Process in parallel
    let products: Vec<PcrProduct> = records
        .par_iter()
        .flat_map(|record| {
            let result: Vec<_> = primers.iter().flat_map(|primer| {
                find_pcr_products(record, primer, config)
            }).collect();
            
            // Update progress
            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(ref pb) = pb {
                pb.set_position(done as u64);
            }
            
            result
        })
        .collect();

    // Finish progress bar
    if let Some(pb) = pb {
        pb.finish_with_message("Search complete");
    }

    // Re-assign case numbers globally
    let mut final_products = products;
    for (i, product) in final_products.iter_mut().enumerate() {
        product.case_number = i + 1;
    }

    final_products
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

        assert_eq!(product.header(), "chr1_amplicon:16S:1");
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

        assert_eq!(product.header(), "chr1_amplicon:16S_rc:2");
    }
}
