use anyhow::{Context, Result, bail};
use std::path::Path;

/// Represents a primer pair for PCR
#[derive(Debug, Clone)]
pub struct PrimerPair {
    /// Name/identifier for this primer pair
    pub name: String,
    /// Forward primer sequence (5' to 3')
    pub forward: Vec<u8>,
    /// Reverse primer sequence (5' to 3')
    pub reverse: Vec<u8>,
}

impl PrimerPair {
    /// Create a new primer pair
    pub fn new(
        name: impl Into<String>,
        forward: impl AsRef<[u8]>,
        reverse: impl AsRef<[u8]>,
    ) -> Result<Self> {
        let name = name.into();
        let forward = forward.as_ref().to_ascii_uppercase();
        let reverse = reverse.as_ref().to_ascii_uppercase();

        // Validate sequences contain only valid IUPAC characters
        validate_iupac_sequence(&forward, &name, "forward")?;
        validate_iupac_sequence(&reverse, &name, "reverse")?;

        Ok(Self {
            name,
            forward,
            reverse,
        })
    }

    /// Get forward primer length
    #[must_use] 
    pub fn forward_len(&self) -> usize {
        self.forward.len()
    }

    /// Get reverse primer length
    #[must_use] 
    pub fn reverse_len(&self) -> usize {
        self.reverse.len()
    }
}

/// Valid IUPAC nucleotide codes
const IUPAC_CODES: &[u8] = b"ACGTRYSWKMBDHVN";

/// Validate that a sequence contains only valid IUPAC characters
fn validate_iupac_sequence(seq: &[u8], primer_name: &str, direction: &str) -> Result<()> {
    for (i, &base) in seq.iter().enumerate() {
        if !IUPAC_CODES.contains(&base) {
            bail!(
                "Invalid character '{}' at position {} in {} primer '{}'. \
                 Valid IUPAC codes are: A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N",
                base as char,
                i + 1,
                direction,
                primer_name
            );
        }
    }
    Ok(())
}

/// Parse primers from either a CLI argument string or a CSV file
pub fn parse_primers(input: &str) -> Result<Vec<PrimerPair>> {
    let path = Path::new(input);

    // Check if input is a file path (CSV)
    if path.exists() && path.is_file() {
        parse_primers_from_csv(path)
    } else {
        // Try to parse as inline primer specification
        parse_primers_from_string(input)
    }
}

/// Parse primers from a CSV file
/// Expected format: name,forward,reverse (with header)
fn parse_primers_from_csv(path: &Path) -> Result<Vec<PrimerPair>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)
        .with_context(|| format!("Failed to open primer CSV file: {}", path.display()))?;

    let mut primers = Vec::new();

    for (i, result) in reader.records().enumerate() {
        let record = result.with_context(|| format!("Error reading CSV row {}", i + 2))?;

        if record.len() < 3 {
            bail!(
                "CSV row {} has {} columns, expected at least 3 (name, forward, reverse)",
                i + 2,
                record.len()
            );
        }

        let name = record.get(0).unwrap_or("").trim();
        let forward = record.get(1).unwrap_or("").trim();
        let reverse = record.get(2).unwrap_or("").trim();

        if name.is_empty() || forward.is_empty() || reverse.is_empty() {
            bail!("CSV row {} has empty fields", i + 2);
        }

        let primer = PrimerPair::new(name, forward.as_bytes(), reverse.as_bytes())
            .with_context(|| format!("Error parsing primer at CSV row {}", i + 2))?;

        primers.push(primer);
    }

    if primers.is_empty() {
        bail!("No primers found in CSV file: {}", path.display());
    }

    log::info!("Loaded {} primer pair(s) from CSV", primers.len());
    Ok(primers)
}

/// Parse primers from a string specification
/// Format: "name:forward:reverse" or multiple separated by semicolons
/// Example: "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT"
fn parse_primers_from_string(input: &str) -> Result<Vec<PrimerPair>> {
    let mut primers = Vec::new();

    // Split by semicolon for multiple primers
    for (i, spec) in input.split(';').enumerate() {
        let spec = spec.trim();
        if spec.is_empty() {
            continue;
        }

        let parts: Vec<&str> = spec.split(':').collect();

        if parts.len() != 3 {
            bail!(
                "Invalid primer specification '{}' (primer {}). \
                 Expected format: 'name:forward:reverse'",
                spec,
                i + 1
            );
        }

        let name = parts[0].trim();
        let forward = parts[1].trim();
        let reverse = parts[2].trim();

        if name.is_empty() || forward.is_empty() || reverse.is_empty() {
            bail!(
                "Primer specification '{spec}' has empty fields. \
                 Expected format: 'name:forward:reverse'"
            );
        }

        let primer = PrimerPair::new(name, forward.as_bytes(), reverse.as_bytes())
            .with_context(|| format!("Error parsing primer specification '{spec}'"))?;

        primers.push(primer);
    }

    if primers.is_empty() {
        bail!(
            "No valid primers found in input '{input}'. \
             Expected format: 'name:forward:reverse' or path to CSV file"
        );
    }

    log::info!("Parsed {} primer pair(s) from command line", primers.len());
    Ok(primers)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_primer() {
        let primer = PrimerPair::new("test", b"ACGT", b"TGCA").unwrap();
        assert_eq!(primer.name, "test");
        assert_eq!(primer.forward, b"ACGT");
        assert_eq!(primer.reverse, b"TGCA");
    }

    #[test]
    fn test_iupac_primer() {
        let primer = PrimerPair::new("iupac", b"ACGTRYSWKMBDHVN", b"ACGT").unwrap();
        assert_eq!(primer.forward, b"ACGTRYSWKMBDHVN");
    }

    #[test]
    fn test_lowercase_conversion() {
        let primer = PrimerPair::new("test", b"acgt", b"tgca").unwrap();
        assert_eq!(primer.forward, b"ACGT");
        assert_eq!(primer.reverse, b"TGCA");
    }

    #[test]
    fn test_invalid_character() {
        let result = PrimerPair::new("test", b"ACGTX", b"TGCA");
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("Invalid character 'X'"));
    }

    #[test]
    fn test_parse_single_primer_string() {
        let primers = parse_primers_from_string("test:ACGT:TGCA").unwrap();
        assert_eq!(primers.len(), 1);
        assert_eq!(primers[0].name, "test");
    }

    #[test]
    fn test_parse_multiple_primer_string() {
        let primers = parse_primers_from_string("p1:ACGT:TGCA;p2:AAAA:TTTT").unwrap();
        assert_eq!(primers.len(), 2);
        assert_eq!(primers[0].name, "p1");
        assert_eq!(primers[1].name, "p2");
    }

    #[test]
    fn test_invalid_primer_format() {
        let result = parse_primers_from_string("invalid");
        assert!(result.is_err());
    }
}
