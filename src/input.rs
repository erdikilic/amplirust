use anyhow::{Context, Result, bail};
use libdeflater::{DecompressionError, Decompressor};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

/// A sequence record from FASTA input
#[derive(Debug, Clone)]
pub struct SequenceRecord {
    /// Header/ID of the sequence (without >)
    pub header: String,
    /// The sequence data (uppercase)
    pub sequence: Vec<u8>,
    /// Source file path
    pub source_file: PathBuf,
}

/// Expand input patterns to a list of files
/// Supports: single files, comma-separated lists, glob patterns
pub fn expand_input_patterns(patterns: &[String]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for pattern in patterns {
        let pattern = pattern.trim();
        if pattern.is_empty() {
            continue;
        }

        // Check if it's a glob pattern
        if pattern.contains('*') || pattern.contains('?') || pattern.contains('[') {
            let matches: Vec<_> = glob::glob(pattern)
                .with_context(|| format!("Invalid glob pattern: {}", pattern))?
                .collect();

            if matches.is_empty() {
                log::warn!("Glob pattern '{}' matched no files", pattern);
            }

            for entry in matches {
                let path = entry.with_context(|| format!("Error reading glob match for '{}'", pattern))?;
                if path.is_file() {
                    files.push(path);
                }
            }
        } else {
            // It's a regular file path
            let path = PathBuf::from(pattern);
            if !path.exists() {
                bail!("Input file not found: {}", path.display());
            }
            if !path.is_file() {
                bail!("Input path is not a file: {}", path.display());
            }
            files.push(path);
        }
    }

    if files.is_empty() {
        bail!("No input files found");
    }

    // Remove duplicates while preserving order
    let mut seen = std::collections::HashSet::new();
    files.retain(|f| seen.insert(f.clone()));

    log::info!("Found {} input file(s)", files.len());
    Ok(files)
}

/// Check if a file is gzip compressed based on extension
pub fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

/// Read and decompress a gzip file using libdeflater
fn read_gzipped_file(path: &Path) -> Result<Vec<u8>> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path.display()))?;
    
    let mut compressed = Vec::new();
    file.read_to_end(&mut compressed)
        .with_context(|| format!("Failed to read file: {}", path.display()))?;

    // Estimate decompressed size (typically 3-10x for text)
    let estimated_size = compressed.len() * 6;
    let mut decompressed = vec![0u8; estimated_size];
    
    let mut decompressor = Decompressor::new();
    
    // Try decompression, grow buffer if needed
    loop {
        match decompressor.gzip_decompress(&compressed, &mut decompressed) {
            Ok(actual_size) => {
                decompressed.truncate(actual_size);
                return Ok(decompressed);
            }
            Err(DecompressionError::InsufficientSpace) => {
                // Double the buffer and retry
                let new_size = decompressed.len() * 2;
                decompressed.resize(new_size, 0);
            }
            Err(e) => {
                bail!("Failed to decompress {}: {:?}", path.display(), e);
            }
        }
    }
}

/// Read a plain (uncompressed) file
fn read_plain_file(path: &Path) -> Result<Vec<u8>> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path.display()))?;
    
    let mut contents = Vec::new();
    file.read_to_end(&mut contents)
        .with_context(|| format!("Failed to read file: {}", path.display()))?;
    
    Ok(contents)
}

/// Read all sequences from a single FASTA file
pub fn read_sequences_from_file(path: &Path) -> Result<Vec<SequenceRecord>> {
    // Read file contents (decompress if gzipped)
    let contents = if is_gzipped(path) {
        read_gzipped_file(path)?
    } else {
        read_plain_file(path)?
    };

    // Parse FASTA using the selected parser
    parse_fasta(&contents, path)
}

// ============================================================================
// FASTA Parser: seq_io implementation
// ============================================================================
#[cfg(feature = "parser_seqio")]
fn parse_fasta(data: &[u8], source: &Path) -> Result<Vec<SequenceRecord>> {
    use seq_io::fasta::{Reader, Record};
    
    log::trace!("Using seq_io parser");
    let mut reader = Reader::new(data);
    let mut records = Vec::new();

    while let Some(result) = reader.next() {
        let record = result.with_context(|| format!("Error reading FASTA from {}", source.display()))?;
        
        let header = record.id()
            .with_context(|| "Invalid UTF-8 in FASTA header")?
            .to_string();
        
        // Collect sequence lines and uppercase
        let sequence: Vec<u8> = record.seq_lines()
            .flat_map(|line| line.iter().copied())
            .map(|b| b.to_ascii_uppercase())
            .collect();

        records.push(SequenceRecord {
            header,
            sequence,
            source_file: source.to_path_buf(),
        });
    }

    log::debug!("Read {} sequences from {}", records.len(), source.display());
    Ok(records)
}

// ============================================================================
// FASTA Parser: needletail implementation
// ============================================================================
#[cfg(feature = "parser_needletail")]
fn parse_fasta(data: &[u8], source: &Path) -> Result<Vec<SequenceRecord>> {
    use needletail::parse_fastx_reader;
    
    log::trace!("Using needletail parser");
    let mut reader = parse_fastx_reader(data)
        .with_context(|| format!("Error opening FASTA from {}", source.display()))?;
    let mut records = Vec::new();

    while let Some(result) = reader.next() {
        let record = result.with_context(|| format!("Error reading FASTA from {}", source.display()))?;
        
        // Get header (needletail includes '>' so we need to handle it)
        let header = std::str::from_utf8(record.id())
            .with_context(|| "Invalid UTF-8 in FASTA header")?
            .to_string();
        
        // Get sequence and uppercase
        let sequence: Vec<u8> = record.seq()
            .iter()
            .copied()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        records.push(SequenceRecord {
            header,
            sequence,
            source_file: source.to_path_buf(),
        });
    }

    log::debug!("Read {} sequences from {}", records.len(), source.display());
    Ok(records)
}

// ============================================================================
// Compile-time check: ensure exactly one parser is selected
// ============================================================================
#[cfg(not(any(feature = "parser_seqio", feature = "parser_needletail")))]
compile_error!("No FASTA parser selected! Enable either 'parser_seqio' or 'parser_needletail' feature.");

#[cfg(all(feature = "parser_seqio", feature = "parser_needletail"))]
compile_error!("Multiple FASTA parsers selected! Enable only one of 'parser_seqio' or 'parser_needletail'.");

/// Return the name of the active FASTA parser
pub fn parser_name() -> &'static str {
    #[cfg(feature = "parser_seqio")]
    {
        "seq_io"
    }
    #[cfg(feature = "parser_needletail")]
    {
        "needletail"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_gzipped() {
        assert!(is_gzipped(Path::new("test.fasta.gz")));
        assert!(is_gzipped(Path::new("test.fa.GZ")));
        assert!(!is_gzipped(Path::new("test.fasta")));
    }

    #[test]
    fn test_parse_fasta() {
        let fasta = b">seq1\nACGT\nTGCA\n>seq2\nAAAA\n";
        let records = parse_fasta(fasta, Path::new("test.fa")).unwrap();
        
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "seq1");
        assert_eq!(records[0].sequence, b"ACGTTGCA");
        assert_eq!(records[1].header, "seq2");
        assert_eq!(records[1].sequence, b"AAAA");
    }

    #[test]
    fn test_parse_fasta_multiline() {
        let fasta = b">long_seq\nACGT\nACGT\nACGT\nACGT\n";
        let records = parse_fasta(fasta, Path::new("test.fa")).unwrap();
        
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, b"ACGTACGTACGTACGT");
    }

    #[test]
    fn test_uppercase_conversion() {
        let fasta = b">seq\nacgt\n";
        let records = parse_fasta(fasta, Path::new("test.fa")).unwrap();
        
        assert_eq!(records[0].sequence, b"ACGT");
    }

    #[test]
    fn test_parser_name() {
        let name = parser_name();
        assert!(name == "seq_io" || name == "needletail");
    }
}
