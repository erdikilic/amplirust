use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use gzp::deflate::Bgzf;
use gzp::par::decompress::ParDecompressBuilder;
use std::fs::File;
use std::io::{BufReader, Read};
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
                let path =
                    entry.with_context(|| format!("Error reading glob match for '{}'", pattern))?;
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

/// Check if a gzip file is in BGZF format by examining header.
/// BGZF files have a specific subfield in the gzip header:
/// - Bytes 12-13: SI1=0x42 ('B'), SI2=0x43 ('C')
/// This "BC" field contains the compressed block size.
fn is_bgzf(path: &Path) -> Result<bool> {
    let mut file = File::open(path)?;
    let mut header = [0u8; 18];

    if file.read_exact(&mut header).is_err() {
        return Ok(false);
    }

    // Check gzip magic (1f 8b) and extra field flag (bit 2 of FLG)
    if header[0] != 0x1f || header[1] != 0x8b || (header[3] & 0x04) == 0 {
        return Ok(false);
    }

    // Check for BGZF subfield identifier "BC" at bytes 12-13
    // XLEN is at bytes 10-11 (little-endian), subfield starts at 12
    Ok(header[12] == b'B' && header[13] == b'C')
}

/// Read and decompress a gzip file using flate2's MultiGzDecoder
/// MultiGzDecoder handles concatenated/multi-member gzip files correctly
fn read_gzipped_file(path: &Path) -> Result<Vec<u8>> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;

    let mut decoder = MultiGzDecoder::new(file);
    let mut decompressed = Vec::new();
    decoder
        .read_to_end(&mut decompressed)
        .with_context(|| format!("Failed to decompress {}", path.display()))?;

    Ok(decompressed)
}

/// Read a plain (uncompressed) file
fn read_plain_file(path: &Path) -> Result<Vec<u8>> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;

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
// Streaming FASTA Reader (for memory-efficient single-file processing)
// ============================================================================

/// Create a streaming FASTA reader that yields sequences lazily.
/// This enables producer-consumer parallelism without loading entire file into memory.
/// For BGZF-compressed files, uses parallel decompression via gzp.
#[cfg(feature = "parser_seqio")]
pub fn streaming_fasta_reader(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<SequenceRecord>> + Send>> {
    let source_file = path.to_path_buf();

    if is_gzipped(path) {
        if is_bgzf(path)? {
            // BGZF: use parallel decompression
            log::debug!("Detected BGZF format, using parallel decompression");
            let file = File::open(path)
                .with_context(|| format!("Failed to open file: {}", path.display()))?;
            let decoder = ParDecompressBuilder::<Bgzf>::new().from_reader(file);
            let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
            Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file)))
        } else {
            // Standard gzip: use single-threaded MultiGzDecoder
            log::debug!("Standard gzip format, using single-threaded decompression");
            let file = File::open(path)
                .with_context(|| format!("Failed to open file: {}", path.display()))?;
            let decoder = MultiGzDecoder::new(file);
            let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
            Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file)))
        }
    } else {
        let file =
            File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
        let buf_reader = BufReader::with_capacity(64 * 1024, file);
        Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file)))
    }
}

#[cfg(feature = "parser_seqio")]
struct StreamingFastaIter<R: Read> {
    reader: seq_io::fasta::Reader<BufReader<R>>,
    source_file: PathBuf,
}

#[cfg(feature = "parser_seqio")]
impl<R: Read> StreamingFastaIter<R> {
    fn new(buf_reader: BufReader<R>, source_file: PathBuf) -> Self {
        Self {
            reader: seq_io::fasta::Reader::new(buf_reader),
            source_file,
        }
    }
}

#[cfg(feature = "parser_seqio")]
impl<R: Read + Send> Iterator for StreamingFastaIter<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        use seq_io::fasta::Record;
        match self.reader.next() {
            Some(Ok(record)) => {
                let header = match record.id() {
                    Ok(id) => id.to_string(),
                    Err(e) => return Some(Err(anyhow::anyhow!("Invalid UTF-8 in header: {}", e))),
                };
                let mut sequence = record.full_seq().into_owned();
                sequence.make_ascii_uppercase();

                Some(Ok(SequenceRecord {
                    header,
                    sequence,
                    source_file: self.source_file.clone(),
                }))
            }
            Some(Err(e)) => Some(Err(anyhow::anyhow!("FASTA parse error: {}", e))),
            None => None,
        }
    }
}

/// Streaming FASTA reader for needletail parser
#[cfg(feature = "parser_needletail")]
pub fn streaming_fasta_reader(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<SequenceRecord>> + Send>> {
    // For needletail, we fall back to loading the file since it doesn't support
    // streaming from a generic Read trait as cleanly as seq_io
    let contents = if is_gzipped(path) {
        read_gzipped_file(path)?
    } else {
        read_plain_file(path)?
    };
    let records = parse_fasta(&contents, path)?;
    Ok(Box::new(records.into_iter().map(Ok)))
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
        let record =
            result.with_context(|| format!("Error reading FASTA from {}", source.display()))?;

        let header = record
            .id()
            .with_context(|| "Invalid UTF-8 in FASTA header")?
            .to_string();

        // Use full_seq() which handles multi-line sequences efficiently
        let mut sequence = record.full_seq().into_owned();
        sequence.make_ascii_uppercase();

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
        let record =
            result.with_context(|| format!("Error reading FASTA from {}", source.display()))?;

        // Get header (needletail includes '>' so we need to handle it)
        let header = std::str::from_utf8(record.id())
            .with_context(|| "Invalid UTF-8 in FASTA header")?
            .to_string();

        // Get sequence and uppercase
        let sequence: Vec<u8> = record
            .seq()
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
compile_error!(
    "No FASTA parser selected! Enable either 'parser_seqio' or 'parser_needletail' feature."
);

#[cfg(all(feature = "parser_seqio", feature = "parser_needletail"))]
compile_error!(
    "Multiple FASTA parsers selected! Enable only one of 'parser_seqio' or 'parser_needletail'."
);

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
    fn test_is_bgzf_standard_gzip() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Create a standard gzip file (not BGZF)
        let mut temp = NamedTempFile::new().unwrap();
        let encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        let mut encoder = encoder;
        encoder.write_all(b">seq\nACGT\n").unwrap();
        let compressed = encoder.finish().unwrap();
        temp.write_all(&compressed).unwrap();
        temp.flush().unwrap();

        // Standard gzip should NOT be detected as BGZF
        assert!(!is_bgzf(temp.path()).unwrap());
    }

    #[test]
    fn test_is_bgzf_too_small() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // File too small to be BGZF
        let mut temp = NamedTempFile::new().unwrap();
        temp.write_all(b"tiny").unwrap();
        temp.flush().unwrap();

        assert!(!is_bgzf(temp.path()).unwrap());
    }

    #[test]
    fn test_is_bgzf_plain_file() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Plain FASTA file (not gzip at all)
        let mut temp = NamedTempFile::new().unwrap();
        temp.write_all(b">seq\nACGT\n").unwrap();
        temp.flush().unwrap();

        assert!(!is_bgzf(temp.path()).unwrap());
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

    #[test]
    fn test_nonexistent_file_error() {
        let result = read_sequences_from_file(Path::new("/nonexistent/path/to/file.fasta"));
        assert!(result.is_err(), "Should error on nonexistent file");

        let err = result.unwrap_err();
        let err_str = format!("{:?}", err);
        assert!(
            err_str.contains("Failed to open") || err_str.contains("No such file"),
            "Error should indicate file not found"
        );
    }

    #[test]
    fn test_invalid_fasta_format() {
        // Data without > header - should fail or return empty
        let invalid_data = b"ACGTACGT\nTGCATGCA\n";
        let result = parse_fasta(invalid_data, Path::new("invalid.fa"));

        // Behavior depends on parser - some skip invalid, some error
        // At minimum, it should not panic
        match result {
            Ok(records) => {
                // If it succeeds, it should have no valid records
                // (or parser interprets it differently)
                assert!(
                    records.is_empty() || records.iter().all(|r| r.header.is_empty()),
                    "Invalid FASTA should produce no valid records or error"
                );
            }
            Err(_) => {
                // Error is also acceptable
            }
        }
    }

    #[test]
    fn test_empty_fasta_file() {
        // Valid but empty FASTA
        let empty_data = b"";
        let result = parse_fasta(empty_data, Path::new("empty.fa"));

        // Parser behavior differs: seq_io returns Ok([]), needletail returns Err
        // Both are acceptable behaviors for an empty file
        match result {
            Ok(records) => assert!(records.is_empty(), "Empty file should produce no records"),
            Err(_) => {
                // needletail errors on empty files - this is acceptable
            }
        }
    }

    #[test]
    fn test_fasta_with_empty_sequence() {
        // FASTA with header but no sequence
        let fasta = b">empty_seq\n>next_seq\nACGT\n";
        let result = parse_fasta(fasta, Path::new("test.fa")).unwrap();

        // Should have 2 records (behavior depends on parser)
        // At least the second one should have sequence
        let non_empty: Vec<_> = result.iter().filter(|r| !r.sequence.is_empty()).collect();
        assert!(
            !non_empty.is_empty(),
            "Should have at least one non-empty sequence"
        );
    }

    #[test]
    fn test_expand_nonexistent_file_error() {
        let result = expand_input_patterns(&["/nonexistent/file.fasta".to_string()]);
        assert!(result.is_err(), "Should error on nonexistent file");
    }

    #[test]
    fn test_expand_empty_patterns() {
        // Empty patterns should return error
        let result = expand_input_patterns(&[]);
        assert!(result.is_err(), "Should error on no input files");

        let result2 = expand_input_patterns(&["".to_string(), "  ".to_string()]);
        assert!(result2.is_err(), "Should error on empty patterns");
    }

    #[test]
    fn test_glob_deduplication() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test.fasta");
        std::fs::write(&file_path, b">seq\nACGT\n").unwrap();

        // Same file referenced twice should be deduplicated
        let patterns = vec![
            file_path.to_string_lossy().to_string(),
            file_path.to_string_lossy().to_string(),
        ];

        let files = expand_input_patterns(&patterns).unwrap();
        assert_eq!(files.len(), 1, "Duplicate files should be deduplicated");
    }
}
