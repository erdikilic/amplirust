use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use gzp::deflate::Bgzf;
use gzp::par::decompress::ParDecompressBuilder;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use crate::errors::ValidationError;

/// Uppercase a sequence in-place, but only if it contains lowercase bytes.
/// Most NCBI sequences are already uppercase, so the scan is a fast early exit.
#[inline]
fn ensure_uppercase(seq: &mut [u8]) {
    if seq.iter().any(u8::is_ascii_lowercase) {
        seq.make_ascii_uppercase();
    }
}

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

/// Supported input file formats
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InputFormat {
    Fasta,
    Genbank,
}

/// Detect the input format from file extension.
/// Returns `None` for unrecognized extensions.
/// Strips `.gz` suffix before checking the inner extension.
#[must_use]
pub fn detect_format(path: &Path) -> Option<InputFormat> {
    // Strip .gz if present to get to the real extension
    let path = if is_gzipped(path) {
        Path::new(path.file_stem()?)
    } else {
        path
    };

    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    match ext.as_str() {
        "fasta" | "fa" | "fna" | "ffn" | "faa" | "frn" | "fas" => Some(InputFormat::Fasta),
        "gb" | "gbk" | "gbff" | "genbank" | "gbf" => Some(InputFormat::Genbank),
        _ => None,
    }
}

/// Expand input patterns to a list of files.
///
/// Supports: single files, comma-separated lists, glob patterns.
///
/// # Errors
///
/// Returns an error if no valid input files are found, or if a file path
/// does not exist.
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
                .with_context(|| format!("Invalid glob pattern: {pattern}"))?
                .collect();

            if matches.is_empty() {
                log::warn!("Glob pattern '{pattern}' matched no files");
            }

            for entry in matches {
                let path =
                    entry.with_context(|| format!("Error reading glob match for '{pattern}'"))?;
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

    // Filter out files with unrecognized extensions
    files.retain(|f| {
        if detect_format(f).is_some() {
            true
        } else {
            log::warn!("Skipping '{}': unrecognized file extension", f.display());
            false
        }
    });

    if files.is_empty() {
        bail!("No supported sequence files found");
    }

    log::info!("Found {} input file(s)", files.len());
    Ok(files)
}

/// Check if a file is gzip compressed based on extension
#[must_use]
pub fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

/// Check if a gzip file is in BGZF format by examining header.
/// BGZF files have a specific subfield in the gzip header:
/// - Bytes 12-13: SI1=0x42 ('B'), SI2=0x43 ('C')
///   This "BC" field contains the compressed block size.
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

/// Read and decompress a gzip file using flate2's `MultiGzDecoder`.
///
/// `MultiGzDecoder` handles concatenated/multi-member gzip files correctly.
/// Pre-allocates based on compressed size x estimated compression ratio
/// to avoid repeated reallocation during decompression.
///
/// When `max_size > 0`, wraps the decoder with `Read::take(max_size + 1)` and
/// checks the result: if more than `max_size` bytes were read, the file exceeds
/// the limit. The `+1` trick distinguishes "exactly at limit" (OK) from "exceeds
/// limit" (reject).
fn read_gzipped_file(path: &Path, max_size: u64) -> Result<Vec<u8>> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let compressed_size = file.metadata().map(|m| m.len()).unwrap_or(0) as usize;

    let decoder = MultiGzDecoder::new(file);
    // Adaptive estimation: GenBank metadata compresses differently than raw sequence.
    // Strip .gz to inspect the inner extension for format-aware ratio selection.
    let inner_path = path.file_stem().map(std::path::Path::new);
    let ratio = if inner_path
        .and_then(|p| p.extension())
        .is_some_and(|e| {
            let lower = e.to_ascii_lowercase();
            matches!(lower.to_str(), Some("gb" | "gbk" | "gbff" | "genbank" | "gbf"))
        }) {
        3 // GenBank files: metadata-heavy, lower compression ratio
    } else {
        4 // FASTA/other: raw sequence, higher compression ratio
    };
    let estimated = compressed_size.saturating_mul(ratio);
    let mut decompressed = Vec::with_capacity(estimated);

    if max_size > 0 {
        let mut limited = decoder.take(max_size + 1);
        limited
            .read_to_end(&mut decompressed)
            .with_context(|| format!("Failed to decompress {}", path.display()))?;
        if decompressed.len() as u64 > max_size {
            return Err(ValidationError::DecompressionLimitExceeded {
                path: path.to_path_buf(),
                limit: max_size,
            }
            .into());
        }
    } else {
        let mut decoder = decoder;
        decoder
            .read_to_end(&mut decompressed)
            .with_context(|| format!("Failed to decompress {}", path.display()))?;
    }

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

/// Read all sequences from a single file (FASTA or `GenBank`).
///
/// # Errors
///
/// Returns an error if the file format is unsupported, the file cannot be
/// read, or the contents cannot be parsed.
pub fn read_sequences_from_file(
    path: &Path,
    max_decompression_size: u64,
) -> Result<Vec<SequenceRecord>> {
    let format = detect_format(path)
        .ok_or_else(|| anyhow::anyhow!("Unsupported format: {}", path.display()))?;

    // Read file contents (decompress if gzipped)
    let contents = if is_gzipped(path) {
        read_gzipped_file(path, max_decompression_size)?
    } else {
        read_plain_file(path)?
    };

    match format {
        InputFormat::Fasta => parse_fasta(&contents, path),
        InputFormat::Genbank => Ok(parse_genbank(&contents, path)),
    }
}

// ============================================================================
// GenBank Parser
// ============================================================================

/// Parse `GenBank` records from raw data (uses slice-based fast path)
fn parse_genbank(data: &[u8], source: &Path) -> Vec<SequenceRecord> {
    use crate::genbank::parse_genbank_slice;

    log::trace!("Using slice-based GenBank parser (fast path)");
    let gb_records = parse_genbank_slice(data);
    let mut records = Vec::new();
    let mut unknown_counter = 0usize;

    for rec in gb_records {
        if rec.seq.is_empty() {
            log::warn!(
                "Skipping record with empty sequence in {}",
                source.display()
            );
            continue;
        }

        let header = rec
            .name
            .as_deref()
            .filter(|s| !s.is_empty())
            .or(rec.accession.as_deref().filter(|s| !s.is_empty()))
            .or(rec.definition.as_deref().filter(|s| !s.is_empty()))
            .map_or_else(
                || {
                    unknown_counter += 1;
                    format!("unknown_{unknown_counter}")
                },
                ToString::to_string,
            );

        if rec.is_circular {
            log::info!(
                "Record '{}' in {} has circular topology; use --circular flag for wrap-around products",
                header,
                source.display()
            );
        }

        let mut sequence = rec.seq;
        ensure_uppercase(&mut sequence);

        records.push(SequenceRecord {
            header,
            sequence,
            source_file: source.to_path_buf(),
        });
    }

    log::debug!("Read {} sequences from {}", records.len(), source.display());
    records
}

// ============================================================================
// Streaming FASTA Reader (for memory-efficient single-file processing)
// ============================================================================

/// Create a streaming reader that yields sequences lazily from FASTA or `GenBank` files.
///
/// This enables producer-consumer parallelism without loading entire file into memory.
/// For BGZF-compressed files, uses parallel decompression via gzp.
///
/// Decompression limit not applied to streaming reader: `seq_io` processes
/// one record at a time from the decompression stream, so memory usage is
/// naturally bounded by a single record's size. The limit exists to prevent
/// `read_to_end()` (used by `read_gzipped_file`) from loading an entire
/// decompression bomb into memory, which cannot happen in streaming mode.
///
/// # Errors
///
/// Returns an error if the file format is unsupported or the file cannot be opened.
#[cfg(feature = "parser_seqio")]
pub fn streaming_reader(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<SequenceRecord>> + Send>> {
    let format = detect_format(path)
        .ok_or_else(|| anyhow::anyhow!("Unsupported format: {}", path.display()))?;
    let source_file = path.to_path_buf();

    if is_gzipped(path) {
        if is_bgzf(path)? {
            // BGZF: use parallel decompression
            log::debug!("Detected BGZF format, using parallel decompression");
            let file = File::open(path)
                .with_context(|| format!("Failed to open file: {}", path.display()))?;
            let decoder = ParDecompressBuilder::<Bgzf>::new().from_reader(file);
            let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
            match format {
                InputFormat::Fasta => {
                    Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file)))
                }
                InputFormat::Genbank => {
                    Ok(Box::new(StreamingGenbankIter::new(buf_reader, source_file)))
                }
            }
        } else {
            // Standard gzip: use single-threaded MultiGzDecoder
            log::debug!("Standard gzip format, using single-threaded decompression");
            let file = File::open(path)
                .with_context(|| format!("Failed to open file: {}", path.display()))?;
            let decoder = MultiGzDecoder::new(file);
            let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
            match format {
                InputFormat::Fasta => {
                    Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file)))
                }
                InputFormat::Genbank => {
                    Ok(Box::new(StreamingGenbankIter::new(buf_reader, source_file)))
                }
            }
        }
    } else {
        let file =
            File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
        let buf_reader = BufReader::with_capacity(64 * 1024, file);
        match format {
            InputFormat::Fasta => Ok(Box::new(StreamingFastaIter::new(buf_reader, source_file))),
            InputFormat::Genbank => {
                Ok(Box::new(StreamingGenbankIter::new(buf_reader, source_file)))
            }
        }
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
                    Err(e) => return Some(Err(anyhow::anyhow!("Invalid UTF-8 in header: {e}"))),
                };
                let mut sequence = record.full_seq().into_owned();
                ensure_uppercase(&mut sequence);

                Some(Ok(SequenceRecord {
                    header,
                    sequence,
                    source_file: self.source_file.clone(),
                }))
            }
            Some(Err(e)) => Some(Err(anyhow::anyhow!("FASTA parse error: {e}"))),
            None => None,
        }
    }
}

// ============================================================================
// Streaming GenBank Reader
// ============================================================================

#[cfg(feature = "parser_seqio")]
struct StreamingGenbankIter<R: Read> {
    reader: crate::genbank::GenbankReader<BufReader<R>>,
    source_file: PathBuf,
    unknown_counter: usize,
}

#[cfg(feature = "parser_seqio")]
impl<R: Read> StreamingGenbankIter<R> {
    fn new(buf_reader: BufReader<R>, source_file: PathBuf) -> Self {
        Self {
            reader: crate::genbank::GenbankReader::from_bufreader(buf_reader),
            source_file,
            unknown_counter: 0,
        }
    }
}

#[cfg(feature = "parser_seqio")]
impl<R: Read + Send> Iterator for StreamingGenbankIter<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.reader.next() {
                Some(Ok(rec)) => {
                    if rec.seq.is_empty() {
                        log::warn!(
                            "Skipping record with empty sequence in {}",
                            self.source_file.display()
                        );
                        continue;
                    }

                    let header = rec
                        .name
                        .as_deref()
                        .filter(|s| !s.is_empty())
                        .or(rec.accession.as_deref().filter(|s| !s.is_empty()))
                        .or(rec.definition.as_deref().filter(|s| !s.is_empty()))
                        .map_or_else(
                            || {
                                self.unknown_counter += 1;
                                format!("unknown_{}", self.unknown_counter)
                            },
                            ToString::to_string,
                        );

                    if rec.is_circular {
                        log::info!(
                            "Record '{header}' has circular topology; use --circular flag for wrap-around products"
                        );
                    }

                    let mut sequence = rec.seq;
                    ensure_uppercase(&mut sequence);

                    return Some(Ok(SequenceRecord {
                        header,
                        sequence,
                        source_file: self.source_file.clone(),
                    }));
                }
                Some(Err(e)) => {
                    return Some(Err(anyhow::anyhow!("GenBank parse error: {e}")));
                }
                None => return None,
            }
        }
    }
}

/// Streaming reader for needletail parser (falls back to loading entire file)
#[cfg(feature = "parser_needletail")]
pub fn streaming_reader(
    path: &Path,
    max_decompression_size: u64,
) -> Result<Box<dyn Iterator<Item = Result<SequenceRecord>> + Send>> {
    let format = detect_format(path)
        .ok_or_else(|| anyhow::anyhow!("Unsupported format: {}", path.display()))?;
    let contents = if is_gzipped(path) {
        read_gzipped_file(path, max_decompression_size)?
    } else {
        read_plain_file(path)?
    };
    let records = match format {
        InputFormat::Fasta => parse_fasta(&contents, path)?,
        InputFormat::Genbank => parse_genbank(&contents, path),
    };
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
        ensure_uppercase(&mut sequence);

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

        // Get sequence and lazy uppercase (skip if already uppercase)
        let mut sequence: Vec<u8> = record.seq().to_vec();
        ensure_uppercase(&mut sequence);

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
#[must_use]
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
        let result = read_sequences_from_file(Path::new("/nonexistent/path/to/file.fasta"), 0);
        assert!(result.is_err(), "Should error on nonexistent file");

        let err = result.unwrap_err();
        let err_str = format!("{err:?}");
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
        if let Ok(records) = result {
            // If it succeeds, it should have no valid records
            // (or parser interprets it differently)
            assert!(
                records.is_empty() || records.iter().all(|r| r.header.is_empty()),
                "Invalid FASTA should produce no valid records or error"
            );
        }
        // Error is also acceptable
    }

    #[test]
    fn test_empty_fasta_file() {
        // Valid but empty FASTA
        let empty_data = b"";
        let result = parse_fasta(empty_data, Path::new("empty.fa"));

        // Parser behavior differs: seq_io returns Ok([]), needletail returns Err
        // Both are acceptable behaviors for an empty file
        if let Ok(records) = result {
            assert!(records.is_empty(), "Empty file should produce no records");
        }
        // needletail errors on empty files - this is acceptable
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

        let result2 = expand_input_patterns(&[String::new(), "  ".to_string()]);
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

    // ============================================================================
    // detect_format() tests
    // ============================================================================

    #[test]
    fn test_detect_format_fasta_extensions() {
        for ext in &["fasta", "fa", "fna", "ffn", "faa", "frn", "fas"] {
            let path = PathBuf::from(format!("test.{ext}"));
            assert_eq!(
                detect_format(&path),
                Some(InputFormat::Fasta),
                "Expected FASTA for .{ext}"
            );
        }
    }

    #[test]
    fn test_detect_format_genbank_extensions() {
        for ext in &["gb", "gbk", "gbff", "genbank", "gbf"] {
            let path = PathBuf::from(format!("test.{ext}"));
            assert_eq!(
                detect_format(&path),
                Some(InputFormat::Genbank),
                "Expected GenBank for .{ext}"
            );
        }
    }

    #[test]
    fn test_detect_format_gz_variants() {
        assert_eq!(
            detect_format(Path::new("test.fasta.gz")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            detect_format(Path::new("test.gbk.gz")),
            Some(InputFormat::Genbank)
        );
        assert_eq!(
            detect_format(Path::new("test.gb.gz")),
            Some(InputFormat::Genbank)
        );
        assert_eq!(
            detect_format(Path::new("test.fa.gz")),
            Some(InputFormat::Fasta)
        );
    }

    #[test]
    fn test_detect_format_unknown() {
        assert_eq!(detect_format(Path::new("test.txt")), None);
        assert_eq!(detect_format(Path::new("test.csv")), None);
        assert_eq!(detect_format(Path::new("test.bed")), None);
        assert_eq!(detect_format(Path::new("test.txt.gz")), None);
    }

    #[test]
    fn test_detect_format_case_insensitive() {
        assert_eq!(
            detect_format(Path::new("test.FASTA")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            detect_format(Path::new("test.Gbk")),
            Some(InputFormat::Genbank)
        );
    }

    // ============================================================================
    // expand_input_patterns() filtering tests
    // ============================================================================

    #[test]
    fn test_expand_skips_non_sequence_files() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        std::fs::write(temp_dir.path().join("seq.fasta"), b">s\nACGT\n").unwrap();
        std::fs::write(temp_dir.path().join("notes.txt"), b"notes").unwrap();
        std::fs::write(temp_dir.path().join("data.csv"), b"a,b,c").unwrap();

        let glob_pattern = format!("{}/*", temp_dir.path().display());
        let files = expand_input_patterns(&[glob_pattern]).unwrap();

        assert_eq!(files.len(), 1);
        assert!(files[0].to_string_lossy().contains("seq.fasta"));
    }

    #[test]
    fn test_expand_all_unrecognized_errors() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        std::fs::write(temp_dir.path().join("notes.txt"), b"notes").unwrap();

        let glob_pattern = format!("{}/*", temp_dir.path().display());
        let result = expand_input_patterns(&[glob_pattern]);
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("No supported sequence files found"));
    }

    // ============================================================================
    // GenBank parsing tests
    // ============================================================================

    /// Minimal valid `GenBank` record for testing
    fn minimal_genbank(name: &str, seq: &str) -> String {
        format!(
            "LOCUS       {:<16} {} bp    DNA     linear   UNK \nDEFINITION  Test sequence.\nACCESSION   ACC001\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\n",
            name,
            seq.len(),
            seq.to_lowercase()
        )
    }

    #[test]
    fn test_parse_genbank_single_record() {
        let gb = minimal_genbank("TestSeq", "acgtacgt");
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, "TestSeq");
        assert_eq!(records[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn test_parse_genbank_multi_record() {
        let gb = format!(
            "{}{}",
            minimal_genbank("Seq1", "aaaa"),
            minimal_genbank("Seq2", "cccc")
        );
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "Seq1");
        assert_eq!(records[0].sequence, b"AAAA");
        assert_eq!(records[1].header, "Seq2");
        assert_eq!(records[1].sequence, b"CCCC");
    }

    #[test]
    fn test_parse_genbank_header_fallback_name() {
        // Name is present — use it
        let gb = minimal_genbank("MyLocus", "acgt");
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));
        assert_eq!(records[0].header, "MyLocus");
    }

    #[test]
    fn test_parse_genbank_header_fallback_unknown() {
        // GenBank with empty name field — should fall back through chain
        let gb = "LOCUS                        4 bp    DNA     linear   UNK \nDEFINITION  .\nACCESSION   \nFEATURES             Location/Qualifiers\nORIGIN\n        1 acgt\n//\n";
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));

        // Should use definition "." or fall back to unknown_1
        assert!(!records[0].header.is_empty());
    }

    #[test]
    fn test_parse_genbank_empty_sequence_skipped() {
        // Record with no ORIGIN/sequence data
        let gb = "LOCUS       EmptySeq             0 bp    DNA     linear   UNK \nDEFINITION  Empty.\nACCESSION   E001\nFEATURES             Location/Qualifiers\nORIGIN\n//\n";
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));
        assert_eq!(records.len(), 0, "Empty sequence records should be skipped");
    }

    #[test]
    fn test_parse_genbank_uppercase() {
        let gb = minimal_genbank("Lower", "acgtnnnn");
        let records = parse_genbank(gb.as_bytes(), Path::new("test.gb"));
        assert_eq!(records[0].sequence, b"ACGTNNNN");
    }

    #[test]
    fn test_read_sequences_from_genbank_file() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gb = minimal_genbank("FileTest", "acgtacgt");
        let mut temp = NamedTempFile::with_suffix(".gb").unwrap();
        temp.write_all(gb.as_bytes()).unwrap();
        temp.flush().unwrap();

        let records = read_sequences_from_file(temp.path(), 0).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, "FileTest");
        assert_eq!(records[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn test_streaming_genbank_reader() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gb = format!(
            "{}{}",
            minimal_genbank("Stream1", "aaaa"),
            minimal_genbank("Stream2", "tttt")
        );
        let mut temp = NamedTempFile::with_suffix(".gb").unwrap();
        temp.write_all(gb.as_bytes()).unwrap();
        temp.flush().unwrap();

        let records: Vec<_> = streaming_reader(temp.path())
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "Stream1");
        assert_eq!(records[0].sequence, b"AAAA");
        assert_eq!(records[1].header, "Stream2");
        assert_eq!(records[1].sequence, b"TTTT");
    }
}
