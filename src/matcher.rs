use crate::utils::calculate_identity;
use sassy::{Match, Searcher, Strand, profiles::Iupac};

/// Configuration for primer matching
#[derive(Debug, Clone)]
pub struct MatchConfig {
    /// Maximum edit distance (mismatches + indels)
    pub max_errors: usize,
    /// Minimum identity threshold (0.0 - 1.0)
    pub min_identity: f64,
    /// Whether to search reverse complement strand
    pub search_rc: bool,
}

impl Default for MatchConfig {
    fn default() -> Self {
        Self {
            max_errors: 2,
            min_identity: 0.0,
            search_rc: false,
        }
    }
}

/// A validated primer match with additional metadata
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    /// Start position in the text (0-based)
    pub start: usize,
    /// End position in the text (exclusive)
    pub end: usize,
    /// Edit distance (cost)
    pub edit_distance: usize,
    /// Strand where match was found
    pub strand: Strand,
    /// CIGAR string describing the alignment
    pub cigar: String,
    /// Calculated identity (0.0 - 1.0)
    pub identity: f64,
}

impl PrimerMatch {
    /// Create from a sassy Match, calculating identity
    pub fn from_sassy_match(m: &Match, primer_len: usize) -> Self {
        let cost = m.cost as usize;
        // Use primer length for identity calculation (standard practice)
        let identity = calculate_identity(cost, primer_len);

        Self {
            start: m.text_start,
            end: m.text_end,
            edit_distance: cost,
            strand: m.strand,
            cigar: m.cigar.to_string(),
            identity,
        }
    }

    /// Length of the matched region in the text
    pub fn match_len(&self) -> usize {
        self.end - self.start
    }
}

/// Primer matcher using sassy for approximate string matching
pub struct PrimerMatcher {
    /// Searcher configured for forward-only or forward+RC
    searcher: Searcher<Iupac>,
    /// Configuration
    config: MatchConfig,
}

impl PrimerMatcher {
    /// Create a new primer matcher with the given configuration
    pub fn new(config: MatchConfig) -> Self {
        let searcher = if config.search_rc {
            Searcher::<Iupac>::new_rc()
        } else {
            Searcher::<Iupac>::new_fwd()
        };

        Self { searcher, config }
    }

    /// Find all matches of a primer in a sequence
    /// Returns matches that pass the identity filter
    pub fn find_matches(&mut self, primer: &[u8], sequence: &[u8]) -> Vec<PrimerMatch> {
        if primer.is_empty() || sequence.is_empty() {
            return Vec::new();
        }

        // Use sassy to find all matches within max_errors
        let matches = self
            .searcher
            .search(primer, sequence, self.config.max_errors);

        // Convert and filter by identity
        matches
            .iter()
            .map(|m| PrimerMatch::from_sassy_match(m, primer.len()))
            .filter(|m| m.identity >= self.config.min_identity)
            .collect()
    }

    /// Find all matches, returning raw sassy matches (for debugging)
    pub fn find_matches_raw(&mut self, primer: &[u8], sequence: &[u8]) -> Vec<Match> {
        if primer.is_empty() || sequence.is_empty() {
            return Vec::new();
        }
        self.searcher
            .search(primer, sequence, self.config.max_errors)
    }

    /// Check if search includes reverse complement
    pub fn searches_rc(&self) -> bool {
        self.config.search_rc
    }
}

/// Find forward primer matches in a sequence
pub fn find_forward_matches(
    primer: &[u8],
    sequence: &[u8],
    config: &MatchConfig,
) -> Vec<PrimerMatch> {
    let mut matcher = PrimerMatcher::new(config.clone());
    matcher.find_matches(primer, sequence)
}

/// Find reverse primer matches in a sequence
/// Note: The reverse primer should be searched on the reverse complement
/// or using sassy's RC mode which handles this automatically
pub fn find_reverse_matches(
    primer: &[u8],
    sequence: &[u8],
    config: &MatchConfig,
) -> Vec<PrimerMatch> {
    let mut matcher = PrimerMatcher::new(config.clone());
    matcher.find_matches(primer, sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let config = MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        let matches = matcher.find_matches(b"ACGT", b"AAAACGTAAA");
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].start, 3);
        assert_eq!(matches[0].end, 7);
        assert_eq!(matches[0].edit_distance, 0);
    }

    #[test]
    fn test_mismatch() {
        let config = MatchConfig {
            max_errors: 1,
            min_identity: 0.5,
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        // ACGT vs ACTT (1 mismatch)
        let matches = matcher.find_matches(b"ACGT", b"AAAACTTAAA");
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].edit_distance, 1);
    }

    #[test]
    fn test_rc_mode() {
        let config = MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: true,
        };
        let mut matcher = PrimerMatcher::new(config);

        // ACGT should match both forward and RC (ACGT is palindromic RC)
        // On sequence with ACGT, should find match
        let matches = matcher.find_matches(b"AAAA", b"AAAATTTT");
        // AAAA forward matches at 0, AAAA RC (TTTT) should match at 4
        assert!(!matches.is_empty());
    }

    #[test]
    fn test_identity_filter() {
        let config = MatchConfig {
            max_errors: 2,
            min_identity: 0.9, // Require 90% identity
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        // With 2 mismatches in 10bp = 80% identity, should be filtered
        let matches = matcher.find_matches(b"AAAAAAAAAA", b"AAAAXXAAAA");
        // Depending on sassy's behavior, this might or might not match
        // The filter should remove low-identity matches
        for m in &matches {
            assert!(m.identity >= 0.9);
        }
    }

    #[test]
    fn test_iupac_matching() {
        let config = MatchConfig {
            max_errors: 0,
            min_identity: 0.8,
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        // R matches A or G
        let matches = matcher.find_matches(b"ACRT", b"AAAACGTAAA");
        // Should match ACGT since R matches G
        assert!(!matches.is_empty());
        // assert!(!matches.is_empty() || true); // IUPAC matching behavior depends on sassy
    }

    #[test]
    fn test_identity_exactly_at_threshold() {
        // Test with identity exactly at threshold
        // 8bp primer with 2 mismatches = 6/8 = 0.75 identity
        let config = MatchConfig {
            max_errors: 2,
            min_identity: 0.75, // Exactly 75%
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        // Primer: ACGTACGT (8bp)
        // Sequence has ACTTACTT (2 mismatches: G->T twice)
        // Identity = (8-2)/8 = 0.75
        let matches = matcher.find_matches(b"ACGTACGT", b"AAAACTTACTTAAA");

        // Should include match at exactly 75% identity (filter uses >=)
        assert!(
            !matches.is_empty(),
            "Should include match at exactly threshold identity"
        );
        for m in &matches {
            assert!(
                m.identity >= 0.75 - 0.001,
                "Identity should be at or above threshold"
            );
        }
    }

    #[test]
    fn test_identity_just_below_threshold() {
        // Test with identity just below threshold - should be filtered
        // 8bp primer with 3 mismatches = 5/8 = 0.625 identity
        let config = MatchConfig {
            max_errors: 3,
            min_identity: 0.75, // Require 75%
            search_rc: false,
        };
        let mut matcher = PrimerMatcher::new(config);

        // Primer: ACGTACGT (8bp)
        // Sequence has ATTTACTT (3 mismatches)
        // Identity = (8-3)/8 = 0.625
        let matches = matcher.find_matches(b"ACGTACGT", b"AAAATTTACTTTAAA");

        // All matches should have identity >= 0.75
        for m in &matches {
            assert!(
                m.identity >= 0.75,
                "Should filter matches below threshold"
            );
        }
    }
}
