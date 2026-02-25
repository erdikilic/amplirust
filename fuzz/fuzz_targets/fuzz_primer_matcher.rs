#![no_main]

use libfuzzer_sys::fuzz_target;
use arbitrary::{Arbitrary, Unstructured, Result as ArbResult};

use amplirust::{MatchConfig, PrimerMatcher};

/// Valid IUPAC nucleotide alphabet for primers
const IUPAC_ALPHABET: &[u8] = b"ACGTRYSWKMBDHVN";
/// Standard DNA alphabet for target sequences
const DNA_ALPHABET: &[u8] = b"ACGT";

/// A fuzz input representing a primer matching scenario
#[derive(Debug)]
struct MatcherInput {
    primer: Vec<u8>,
    sequence: Vec<u8>,
    max_errors: usize,
}

impl<'a> Arbitrary<'a> for MatcherInput {
    fn arbitrary(u: &mut Unstructured<'a>) -> ArbResult<Self> {
        // Primer: 5-50 bases from IUPAC alphabet
        let primer_len = u.int_in_range(5..=50)?;
        let primer: Vec<u8> = (0..primer_len)
            .map(|_| {
                let idx = u.int_in_range(0..=(IUPAC_ALPHABET.len() - 1))?;
                Ok(IUPAC_ALPHABET[idx])
            })
            .collect::<ArbResult<_>>()?;

        // Target sequence: 10-500 bases from DNA alphabet
        let seq_len = u.int_in_range(10..=500)?;
        let sequence: Vec<u8> = (0..seq_len)
            .map(|_| {
                let idx = u.int_in_range(0..=(DNA_ALPHABET.len() - 1))?;
                Ok(DNA_ALPHABET[idx])
            })
            .collect::<ArbResult<_>>()?;

        // Max errors: 0-5
        let max_errors = u.int_in_range(0..=5)?;

        Ok(MatcherInput {
            primer,
            sequence,
            max_errors,
        })
    }
}

fuzz_target!(|input: MatcherInput| {
    let config = MatchConfig {
        max_errors: input.max_errors,
        min_identity: 0.0,
        search_rc: true,
    };
    let mut matcher = PrimerMatcher::new(config);
    let matches = matcher.find_matches(&input.primer, &input.sequence);
    for m in &matches {
        let _ = m.start;
        let _ = m.end;
        let _ = m.edit_distance;
        let _ = m.identity;
        let _ = &m.cigar;
    }
});
