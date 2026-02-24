//! Primer pair and pool parsing and validation.
//!
//! Reads primer pairs (or individual pool primers) from CSV files or inline
//! colon-separated strings, validates IUPAC nucleotide characters, and warns
//! about unusual lengths.

use anyhow::{Context, Result, bail};
use std::path::Path;

use crate::errors::ValidationError;

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
    /// Create a new primer pair.
    ///
    /// # Errors
    ///
    /// Returns an error if either primer contains invalid IUPAC characters.
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

/// Represents a single primer for pool mode (all-vs-all matching)
#[derive(Debug, Clone)]
pub struct Primer {
    /// Name/identifier for this primer
    pub name: String,
    /// Primer sequence (5' to 3')
    pub sequence: Vec<u8>,
}

impl Primer {
    /// Create a new pool primer.
    ///
    /// # Errors
    ///
    /// Returns an error if the sequence contains invalid IUPAC characters.
    pub fn new(name: impl Into<String>, sequence: impl AsRef<[u8]>) -> Result<Self> {
        let name = name.into();
        let sequence = sequence.as_ref().to_ascii_uppercase();
        validate_iupac_sequence(&sequence, &name, "pool")?;
        Ok(Self { name, sequence })
    }

    /// Get primer sequence length
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the primer sequence is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Minimum recommended primer length (bp). Primers shorter than this are
/// unusually short for typical PCR and may produce non-specific amplification.
const MIN_PRIMER_LEN: usize = 10;

/// Maximum recommended primer length (bp). Primers longer than this are
/// uncommon and may indicate a data entry error.
const MAX_PRIMER_LEN: usize = 50;

/// Emit `log::warn!` for each primer arm whose length falls outside the
/// recommended range (`MIN_PRIMER_LEN..=MAX_PRIMER_LEN`). This is advisory
/// only -- bioinformatics convention: unusual inputs should not block execution.
pub fn warn_primer_length(primer: &PrimerPair) {
    for (direction, len) in [
        ("forward", primer.forward_len()),
        ("reverse", primer.reverse_len()),
    ] {
        if !(MIN_PRIMER_LEN..=MAX_PRIMER_LEN).contains(&len) {
            log::warn!(
                "Primer '{}' {} arm has unusual length ({} bp); recommended range is {}-{} bp",
                primer.name,
                direction,
                len,
                MIN_PRIMER_LEN,
                MAX_PRIMER_LEN,
            );
        }
    }
}

/// Emit `log::warn!` for a pool primer whose length falls outside the
/// recommended range (`MIN_PRIMER_LEN..=MAX_PRIMER_LEN`).
pub fn warn_pool_primer_length(primer: &Primer) {
    let len = primer.len();
    if !(MIN_PRIMER_LEN..=MAX_PRIMER_LEN).contains(&len) {
        log::warn!(
            "Pool primer '{}' has unusual length ({} bp); recommended range is {}-{} bp",
            primer.name,
            len,
            MIN_PRIMER_LEN,
            MAX_PRIMER_LEN,
        );
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

/// Parse primers from either a CLI argument string or a CSV file.
///
/// # Errors
///
/// Returns an error if the input cannot be parsed as a valid primer specification
/// or if the CSV file cannot be read.
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

/// Parse pool primers from either a CLI argument string or a CSV file.
///
/// Pool format uses 2 fields (name + sequence) instead of 3 (name + forward + reverse).
///
/// # Errors
///
/// Returns an error if the input cannot be parsed as valid pool primers
/// or if the CSV file cannot be read.
pub fn parse_pool_primers(input: &str) -> Result<Vec<Primer>> {
    let path = Path::new(input);

    if path.exists() && path.is_file() {
        parse_pool_primers_from_csv(path)
    } else {
        parse_pool_primers_from_string(input)
    }
}

/// Parse pool primers from a CSV file
/// Expected format: name,sequence (with header, 2 columns)
fn parse_pool_primers_from_csv(path: &Path) -> Result<Vec<Primer>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)
        .with_context(|| format!("Failed to open pool primer CSV file: {}", path.display()))?;

    let headers = reader
        .headers()
        .with_context(|| {
            format!(
                "Failed to read CSV headers from '{}'; expected format: name,sequence",
                path.display()
            )
        })?
        .clone();

    if headers.len() < 2 {
        return Err(ValidationError::CsvFormat {
            path: path.to_path_buf(),
            detail: format!(
                "header has {} column(s), expected at least 2 (name, sequence). Got: {}",
                headers.len(),
                headers.iter().collect::<Vec<_>>().join(","),
            ),
        }
        .into());
    }

    let mut primers = Vec::new();

    for (i, result) in reader.records().enumerate() {
        let record = result.with_context(|| format!("Error reading CSV row {}", i + 2))?;

        if record.len() < 2 {
            return Err(ValidationError::CsvFormat {
                path: path.to_path_buf(),
                detail: format!(
                    "row {} has {} column(s), expected at least 2 (name, sequence)",
                    i + 2,
                    record.len(),
                ),
            }
            .into());
        }

        let name = record.get(0).unwrap_or("").trim();
        let sequence = record.get(1).unwrap_or("").trim();

        if name.is_empty() || sequence.is_empty() {
            bail!("Pool CSV row {} has empty fields", i + 2);
        }

        let primer = Primer::new(name, sequence.as_bytes())
            .with_context(|| format!("Error parsing pool primer at CSV row {}", i + 2))?;

        primers.push(primer);
    }

    if primers.is_empty() {
        bail!("No pool primers found in CSV file: {}", path.display());
    }

    log::info!("Loaded {} pool primer(s) from CSV", primers.len());
    Ok(primers)
}

/// Parse pool primers from a string specification
/// Format: "name:sequence" or multiple separated by semicolons
/// Example: "p1:AGAGTTTGATCMTGGCTCAG;p2:GWATTACCGCGGCKGCTG"
fn parse_pool_primers_from_string(input: &str) -> Result<Vec<Primer>> {
    let mut primers = Vec::new();

    for (i, spec) in input.split(';').enumerate() {
        let spec = spec.trim();
        if spec.is_empty() {
            continue;
        }

        let parts: Vec<&str> = spec.split(':').collect();

        if parts.len() != 2 {
            bail!(
                "Invalid pool primer specification '{}' (primer {}). \
                 Expected format: 'name:sequence'",
                spec,
                i + 1
            );
        }

        let name = parts[0].trim();
        let sequence = parts[1].trim();

        if name.is_empty() || sequence.is_empty() {
            bail!(
                "Pool primer specification '{spec}' has empty fields. \
                 Expected format: 'name:sequence'"
            );
        }

        let primer = Primer::new(name, sequence.as_bytes())
            .with_context(|| format!("Error parsing pool primer specification '{spec}'"))?;

        primers.push(primer);
    }

    if primers.is_empty() {
        bail!(
            "No valid pool primers found in input '{input}'. \
             Expected format: 'name:sequence' or path to CSV file"
        );
    }

    log::info!("Parsed {} pool primer(s) from command line", primers.len());
    Ok(primers)
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

    // Validate CSV headers before iterating records
    let headers = reader
        .headers()
        .with_context(|| {
            format!(
                "Failed to read CSV headers from '{}'; expected format: name,forward,reverse",
                path.display()
            )
        })?
        .clone();

    if headers.len() < 3 {
        return Err(ValidationError::CsvFormat {
            path: path.to_path_buf(),
            detail: format!(
                "header has {} column(s), expected at least 3 (name, forward, reverse). Got: {}",
                headers.len(),
                headers.iter().collect::<Vec<_>>().join(","),
            ),
        }
        .into());
    }

    let mut primers = Vec::new();

    for (i, result) in reader.records().enumerate() {
        let record = result.with_context(|| format!("Error reading CSV row {}", i + 2))?;

        if record.len() < 3 {
            return Err(ValidationError::CsvFormat {
                path: path.to_path_buf(),
                detail: format!(
                    "row {} has {} column(s), expected at least 3 (name, forward, reverse)",
                    i + 2,
                    record.len(),
                ),
            }
            .into());
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

    #[test]
    fn test_csv_too_few_header_columns() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name,forward").unwrap();
        writeln!(temp, "p1,ACGT").unwrap();
        temp.flush().unwrap();

        let result = parse_primers(temp.path().to_str().unwrap());
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("column"),
            "Expected column count error, got: {err}"
        );
    }

    #[test]
    fn test_csv_row_too_few_columns() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name,forward,reverse").unwrap();
        writeln!(temp, "p1,ACGT").unwrap(); // only 2 columns in data row
        temp.flush().unwrap();

        let result = parse_primers(temp.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_csv_empty_fields() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name,forward,reverse").unwrap();
        writeln!(temp, "p1,,TGCA").unwrap(); // empty forward field
        temp.flush().unwrap();

        let result = parse_primers(temp.path().to_str().unwrap());
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("empty"),
            "Expected empty fields error, got: {err}"
        );
    }

    #[test]
    fn test_warn_primer_length_short() {
        // Primers shorter than MIN_PRIMER_LEN (10) should not panic
        let primer = PrimerPair::new("short", b"ACGT", b"TGCA").unwrap();
        warn_primer_length(&primer); // no panic = pass
    }

    #[test]
    fn test_warn_primer_length_long() {
        // Primers longer than MAX_PRIMER_LEN (50) should not panic
        let long_seq = b"ACGT".repeat(15); // 60bp
        let primer = PrimerPair::new("long", &long_seq, &long_seq).unwrap();
        warn_primer_length(&primer); // no panic = pass
    }

    #[test]
    fn test_warn_primer_length_normal() {
        // Primers within recommended range should not panic
        let primer =
            PrimerPair::new("normal", b"ACGTACGTACGTACGTACGT", b"TGCATGCATGCATGCATGCA").unwrap();
        warn_primer_length(&primer); // no panic = pass
    }

    // ── Pool primer tests ─────────────────────────────────────────────────

    #[test]
    fn test_pool_primer_new() {
        let primer = Primer::new("p1", b"ACGT").unwrap();
        assert_eq!(primer.name, "p1");
        assert_eq!(primer.sequence, b"ACGT");
        assert_eq!(primer.len(), 4);
        assert!(!primer.is_empty());
    }

    #[test]
    fn test_pool_primer_lowercase() {
        let primer = Primer::new("p1", b"acgt").unwrap();
        assert_eq!(primer.sequence, b"ACGT");
    }

    #[test]
    fn test_pool_primer_invalid_char() {
        let result = Primer::new("p1", b"ACGTX");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_pool_primers_inline() {
        let primers = parse_pool_primers("p1:ACGT;p2:TGCA").unwrap();
        assert_eq!(primers.len(), 2);
        assert_eq!(primers[0].name, "p1");
        assert_eq!(primers[0].sequence, b"ACGT");
        assert_eq!(primers[1].name, "p2");
        assert_eq!(primers[1].sequence, b"TGCA");
    }

    #[test]
    fn test_parse_pool_primers_single() {
        let primers = parse_pool_primers("p1:ACGT").unwrap();
        assert_eq!(primers.len(), 1);
    }

    #[test]
    fn test_parse_pool_primers_invalid_format() {
        // 3-field format should fail in pool mode
        let result = parse_pool_primers("p1:ACGT:TGCA");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_pool_primers_empty_fields() {
        let result = parse_pool_primers("p1:");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_pool_primers_csv() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name,sequence").unwrap();
        writeln!(temp, "p1,AGAGTTTGATCMTGGCTCAG").unwrap();
        writeln!(temp, "p2,GWATTACCGCGGCKGCTG").unwrap();
        temp.flush().unwrap();

        let primers = parse_pool_primers(temp.path().to_str().unwrap()).unwrap();
        assert_eq!(primers.len(), 2);
        assert_eq!(primers[0].name, "p1");
        assert_eq!(primers[1].name, "p2");
    }

    #[test]
    fn test_parse_pool_primers_csv_too_few_columns() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name").unwrap();
        writeln!(temp, "p1").unwrap();
        temp.flush().unwrap();

        let result = parse_pool_primers(temp.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_pool_primers_csv_empty_fields() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(temp, "name,sequence").unwrap();
        writeln!(temp, "p1,").unwrap();
        temp.flush().unwrap();

        let result = parse_pool_primers(temp.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_warn_pool_primer_length() {
        let short = Primer::new("short", b"ACGT").unwrap();
        warn_pool_primer_length(&short); // no panic = pass

        let normal = Primer::new("normal", b"ACGTACGTACGTACGTACGT").unwrap();
        warn_pool_primer_length(&normal); // no panic = pass
    }
}
