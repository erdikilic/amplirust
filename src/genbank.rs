//! Lightweight `GenBank` format parser.
//!
//! Only extracts fields needed for PCR analysis: LOCUS name, DEFINITION,
//! ACCESSION, topology (circular/linear), and ORIGIN sequence data.
//! All other sections (FEATURES, REFERENCES, COMMENTS, etc.) are skipped.

use anyhow::{Context, Result};
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use crate::errors::ValidationError;

/// Maximum allowed line length in bytes. Protects against pathological input
/// (e.g., a file with no newlines designed to exhaust memory). Normal `GenBank`
/// lines are under 200 bytes; long DEFINITION or COMMENT lines rarely exceed
/// a few hundred. This limit is a safety net, not a format enforcer.
const MAX_LINE_LENGTH: usize = 100_000; // 100 KB

/// Check that a line does not exceed `MAX_LINE_LENGTH`.
///
/// Uses `trim_end()` length because `BufRead::read_line` appends the trailing
/// newline to the buffer, which should not count against the content limit.
fn check_line_length(line: &str) -> Result<()> {
    let content_len = line.trim_end().len();
    if content_len > MAX_LINE_LENGTH {
        return Err(ValidationError::LineTooLong {
            path: PathBuf::from("<genbank input>"),
            len: content_len,
            limit: MAX_LINE_LENGTH,
        }
        .into());
    }
    Ok(())
}

/// A parsed `GenBank` record containing only the fields needed for PCR analysis.
pub struct GenbankRecord {
    /// LOCUS name (first token of the LOCUS line), used as the sequence identifier.
    pub name: Option<String>,
    /// Primary accession number from the ACCESSION line.
    pub accession: Option<String>,
    /// Full organism/molecule description from the DEFINITION line(s).
    pub definition: Option<String>,
    /// Whether the LOCUS line indicates circular topology.
    pub is_circular: bool,
    /// Raw nucleotide sequence bytes extracted from the ORIGIN section.
    pub seq: Vec<u8>,
}

/// Streaming iterator that yields one `GenbankRecord` at a time from any buffered reader.
pub struct GenbankReader<B: BufRead> {
    reader: B,
    line_buf: String,
    /// True if we've already read a LOCUS line into `line_buf` for the next record.
    has_pending_locus: bool,
}

impl<R: Read> GenbankReader<BufReader<R>> {
    /// Create a new reader, wrapping the given `Read` in a `BufReader`.
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            // Pre-allocate to typical GenBank ORIGIN line width (~80 chars:
            // 9-char line number + 6 groups of 10 bases + spaces) to avoid
            // early reallocations during sequential line reads.
            line_buf: String::with_capacity(80),
            has_pending_locus: false,
        }
    }
}

impl<B: BufRead> GenbankReader<B> {
    /// Wrap an existing buffered reader (avoids double-buffering).
    pub fn from_bufreader(reader: B) -> Self {
        Self {
            reader,
            // Pre-allocate to typical GenBank ORIGIN line width (~80 chars)
            line_buf: String::with_capacity(80),
            has_pending_locus: false,
        }
    }
}

impl<B: BufRead> Iterator for GenbankReader<B> {
    type Item = Result<GenbankRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl<B: BufRead> GenbankReader<B> {
    fn read_record(&mut self) -> Result<Option<GenbankRecord>> {
        // Find the LOCUS line (start of a record)
        if !self.has_pending_locus {
            loop {
                self.line_buf.clear();
                let n = self
                    .reader
                    .read_line(&mut self.line_buf)
                    .context("Failed to read GenBank data")?;
                if n == 0 {
                    return Ok(None); // EOF
                }
                check_line_length(&self.line_buf)?;
                if self.line_buf.starts_with("LOCUS") {
                    break;
                }
            }
        }
        self.has_pending_locus = false;

        let mut record = GenbankRecord {
            name: None,
            accession: None,
            definition: None,
            is_circular: false,
            seq: Vec::new(),
        };

        // Parse LOCUS line
        let (name, is_circular) = parse_locus_fields(self.line_buf.as_bytes());
        record.name = name;
        record.is_circular = is_circular;

        // Read subsequent lines until record terminator
        loop {
            self.line_buf.clear();
            let n = self
                .reader
                .read_line(&mut self.line_buf)
                .context("Failed to read GenBank data")?;
            if n == 0 {
                // EOF mid-record — return what we have
                log::warn!(
                    "Truncated GenBank record '{}': EOF reached before record terminator '//'",
                    record.name.as_deref().unwrap_or("<unnamed>")
                );
                break;
            }
            check_line_length(&self.line_buf)?;

            let line = self.line_buf.trim_end();

            if line == "//" {
                break;
            }

            if line.starts_with("DEFINITION") {
                let text = line.get(12..).unwrap_or("").trim();
                let mut def = text.to_string();
                // Read continuation lines (start with whitespace, not a keyword)
                loop {
                    self.line_buf.clear();
                    let n = self
                        .reader
                        .read_line(&mut self.line_buf)
                        .context("Failed to read GenBank data")?;
                    if n == 0 {
                        break;
                    }
                    check_line_length(&self.line_buf)?;
                    if !self.line_buf.starts_with(' ') || is_keyword_bytes(self.line_buf.as_bytes())
                    {
                        // Not a continuation — this line belongs to the next field
                        break;
                    }
                    def.push(' ');
                    def.push_str(self.line_buf.trim());
                }
                record.definition = Some(def);
                // The line_buf now has the next field line; process it in-place
                // by reprocessing from the top of the loop. We need to check it.
                let reprocess = self.line_buf.clone();
                self.line_buf = reprocess;
                // Check if we need to handle this line
                if self.line_buf.trim_end() == "//" {
                    break;
                }
                if self.line_buf.starts_with("ACCESSION") {
                    let text = self.line_buf.trim_end();
                    if let Some(acc) = text
                        .get(12..)
                        .and_then(|s| s.split_whitespace().next())
                        .filter(|a| !a.is_empty())
                    {
                        record.accession = Some(acc.to_string());
                    }
                } else if self.line_buf.starts_with("ORIGIN") {
                    if read_origin_sequence(&mut self.reader, &mut self.line_buf, &mut record)? {
                        self.has_pending_locus = true;
                    }
                    break;
                } else if self.line_buf.starts_with("LOCUS") {
                    // Next record started without terminator
                    self.has_pending_locus = true;
                    break;
                }
                continue;
            }

            if line.starts_with("ACCESSION") {
                if let Some(acc) = line
                    .get(12..)
                    .and_then(|s| s.split_whitespace().next())
                    .filter(|a| !a.is_empty())
                {
                    record.accession = Some(acc.to_string());
                }
                continue;
            }

            if line.starts_with("ORIGIN") {
                if read_origin_sequence(&mut self.reader, &mut self.line_buf, &mut record)? {
                    self.has_pending_locus = true;
                }
                break;
            }

            if line.starts_with("LOCUS") {
                // Next record started without terminator
                self.has_pending_locus = true;
                break;
            }

            // Skip everything else (FEATURES, REFERENCE, COMMENT, SOURCE, VERSION, etc.)
        }

        Ok(Some(record))
    }
}

/// Extract name and circularity from a LOCUS line (byte slice).
///
/// Format: `LOCUS       name           len bp    mol  topology division`
/// The name field starts at column 12 and extends to the first whitespace.
/// Topology is either "circular" or "linear" somewhere on the line.
///
/// Used by both the streaming and slice-based parser paths.
fn parse_locus_fields(line: &[u8]) -> (Option<String>, bool) {
    let name = if line.len() > 12 {
        let rest = &line[12..];
        let name_end = rest
            .iter()
            .position(|&b| b == b' ' || b == b'\t')
            .unwrap_or(rest.len());
        if name_end > 0 {
            Some(String::from_utf8_lossy(&rest[..name_end]).into_owned())
        } else {
            None
        }
    } else {
        None
    };
    let is_circular = line.windows(8).any(|w| w.eq_ignore_ascii_case(b"circular"));
    (name, is_circular)
}

/// Read sequence data from ORIGIN section until `//`, `LOCUS`, or EOF.
/// Returns `true` if the stop was caused by encountering a new LOCUS line
/// (meaning `line_buf` holds the LOCUS line for the next record).
fn read_origin_sequence<B: BufRead>(
    reader: &mut B,
    line_buf: &mut String,
    record: &mut GenbankRecord,
) -> Result<bool> {
    loop {
        line_buf.clear();
        let n = reader
            .read_line(line_buf)
            .context("Failed to read GenBank sequence data")?;
        if n == 0 {
            log::warn!(
                "Truncated GenBank record: EOF reached in ORIGIN section before record terminator '//'",
            );
            return Ok(false); // EOF
        }
        check_line_length(line_buf)?;
        let trimmed = line_buf.trim();
        if trimmed == "//" {
            return Ok(false);
        }
        if line_buf.starts_with("LOCUS") {
            return Ok(true); // Next record encountered
        }
        // Extract only alphabetic characters (skip line numbers and spaces)
        for &b in line_buf.as_bytes() {
            if b.is_ascii_alphabetic() {
                record.seq.push(b);
            }
        }
    }
}

// ============================================================================
// Slice-based fast path for in-memory GenBank data
// ============================================================================

use memchr::memmem;

/// Parse all `GenBank` records from an in-memory byte slice.
///
/// This is a fast path that operates directly on `&[u8]` without
/// `BufReader`/`read_line` overhead. It uses `memchr::memmem` to find record
/// boundaries (`\nLOCUS`) and jump directly to `\nORIGIN`, skipping FEATURES
/// entirely.
#[must_use]
pub fn parse_genbank_slice(data: &[u8]) -> Vec<GenbankRecord> {
    if data.is_empty() {
        return Vec::new();
    }

    let locus_finder = memmem::Finder::new(b"\nLOCUS");
    let origin_finder = memmem::FinderRev::new(b"\nORIGIN");

    // Collect all record start positions.
    // The first record may start at offset 0 (no preceding newline).
    let mut starts: Vec<usize> = Vec::new();
    if data.starts_with(b"LOCUS") {
        starts.push(0);
    }
    for pos in locus_finder.find_iter(data) {
        starts.push(pos + 1); // skip the '\n'
    }

    if starts.is_empty() {
        return Vec::new();
    }

    let mut records = Vec::with_capacity(starts.len());
    for i in 0..starts.len() {
        let start = starts[i];
        let end = if i + 1 < starts.len() {
            starts[i + 1]
        } else {
            data.len()
        };
        let slice = &data[start..end];
        records.push(parse_single_record_slice(slice, &origin_finder));
    }

    records
}

/// Parse a single `GenBank` record from a byte slice that starts with `LOCUS`.
fn parse_single_record_slice(data: &[u8], origin_finder: &memmem::FinderRev) -> GenbankRecord {
    let mut record = GenbankRecord {
        name: None,
        accession: None,
        definition: None,
        is_circular: false,
        seq: Vec::new(),
    };

    // Find the end of the LOCUS line.
    let locus_end = memchr::memchr(b'\n', data).unwrap_or(data.len());
    let locus_line = &data[..locus_end];
    let (name, is_circular) = parse_locus_fields(locus_line);
    record.name = name;
    record.is_circular = is_circular;

    // Use rfind to locate ORIGIN from the end of the record — this skips over
    // FEATURES (typically the largest section) instead of scanning through it.
    let origin_pos = origin_finder.rfind(data);

    // Parse header fields from just the header portion (before FEATURES or ORIGIN).
    // Use origin_pos as an upper bound to avoid scanning into FEATURES.
    let header_end = find_header_end(data, locus_end, origin_pos);
    if header_end > locus_end + 1 {
        parse_header_fields(&data[locus_end + 1..header_end], &mut record);
    }

    // Extract ORIGIN sequence if present.
    if let Some(opos) = origin_pos {
        // opos points to the '\n' before "ORIGIN"; the ORIGIN line starts at opos+1.
        // Find the end of the ORIGIN header line.
        let origin_line_end =
            opos + 1 + memchr::memchr(b'\n', &data[opos + 1..]).unwrap_or(data[opos + 1..].len());
        let seq_region = &data[origin_line_end..];
        // Find the record terminator "//" within the sequence region.
        let term_pos = memmem::find(seq_region, b"\n//").map_or(seq_region.len(), |p| p + 1);
        extract_origin_fast(&seq_region[..term_pos], &mut record.seq);
    }

    record
}

/// Find where the header section ends (before FEATURES or ORIGIN).
/// If `origin_pos` is known, uses it as an upper bound to avoid scanning into FEATURES.
fn find_header_end(data: &[u8], after_locus: usize, origin_pos: Option<usize>) -> usize {
    // If we know where ORIGIN is, limit the scan region.
    let limit = origin_pos.map_or(data.len(), |p| p + 1);
    // Scan line-by-line from after the LOCUS line.
    let mut pos = after_locus + 1;
    while pos < limit {
        if data[pos..].starts_with(b"FEATURES")
            || data[pos..].starts_with(b"ORIGIN")
            || data[pos..].starts_with(b"//")
        {
            return pos;
        }
        // Advance to next line.
        match memchr::memchr(b'\n', &data[pos..limit]) {
            Some(nl) => pos += nl + 1,
            None => return limit,
        }
    }
    limit
}

/// Parse DEFINITION and ACCESSION from the header portion (bytes between LOCUS and FEATURES/ORIGIN).
fn parse_header_fields(data: &[u8], record: &mut GenbankRecord) {
    let mut pos = 0;
    while pos < data.len() {
        let line_end = pos + memchr::memchr(b'\n', &data[pos..]).unwrap_or(data[pos..].len());
        let line = &data[pos..line_end];

        if line.starts_with(b"DEFINITION") {
            let text_start = if line.len() > 12 { 12 } else { line.len() };
            let mut def = trim_bytes(&line[text_start..]).to_vec();

            // Read continuation lines.
            let mut cpos = line_end + 1;
            while cpos < data.len() {
                let cline_end =
                    cpos + memchr::memchr(b'\n', &data[cpos..]).unwrap_or(data[cpos..].len());
                let cline = &data[cpos..cline_end];
                if cline.is_empty() || !cline[0].is_ascii_whitespace() || is_keyword_bytes(cline) {
                    break;
                }
                def.push(b' ');
                def.extend_from_slice(trim_bytes(cline));
                cpos = cline_end + 1;
            }
            record.definition = Some(String::from_utf8_lossy(&def).into_owned());
            pos = cpos;
            continue;
        }

        if line.starts_with(b"ACCESSION") {
            let text_start = if line.len() > 12 { 12 } else { line.len() };
            let rest = trim_bytes(&line[text_start..]);
            // Take first whitespace-delimited token.
            let acc_end = rest
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(rest.len());
            if acc_end > 0 {
                record.accession = Some(String::from_utf8_lossy(&rest[..acc_end]).into_owned());
            }
            pos = line_end + 1;
            continue;
        }

        pos = line_end + 1;
    }
}

/// Extract sequence data from ORIGIN section lines.
///
/// Pre-allocates capacity and uses `memchr` for fast newline detection.
/// Each ORIGIN line has a ~10 char line-number prefix followed by sequence groups.
fn extract_origin_fast(data: &[u8], seq: &mut Vec<u8>) {
    // Rough estimate: ~80% of bytes are sequence characters.
    seq.reserve(data.len() * 4 / 5);

    let mut pos = 0;
    while pos < data.len() {
        // Skip leading whitespace / newlines.
        if data[pos] == b'\n' || data[pos] == b'\r' {
            pos += 1;
            continue;
        }

        // Find end of this line.
        let line_end = pos + memchr::memchr(b'\n', &data[pos..]).unwrap_or(data[pos..].len());
        let line = &data[pos..line_end];

        // Skip the line-number prefix (first 10 bytes typically).
        // We just filter for alphabetic bytes which is robust regardless of prefix length.
        for &b in line {
            if b.is_ascii_alphabetic() {
                seq.push(b);
            }
        }

        pos = line_end + 1;
    }
}

/// Trim leading and trailing ASCII whitespace from a byte slice.
fn trim_bytes(data: &[u8]) -> &[u8] {
    let start = data
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .unwrap_or(data.len());
    let end = data
        .iter()
        .rposition(|b| !b.is_ascii_whitespace())
        .map_or(start, |p| p + 1);
    &data[start..end]
}

/// Check if a byte line starts with a known `GenBank` top-level keyword.
fn is_keyword_bytes(line: &[u8]) -> bool {
    const KEYWORDS: &[&[u8]] = &[
        b"LOCUS",
        b"DEFINITION",
        b"ACCESSION",
        b"VERSION",
        b"DBLINK",
        b"KEYWORDS",
        b"SOURCE",
        b"REFERENCE",
        b"COMMENT",
        b"FEATURES",
        b"BASE COUNT",
        b"ORIGIN",
        b"CONTIG",
    ];
    KEYWORDS.iter().any(|kw| line.starts_with(kw))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn parse_all(data: &str) -> Vec<GenbankRecord> {
        GenbankReader::new(Cursor::new(data.as_bytes()))
            .collect::<Result<Vec<_>>>()
            .unwrap()
    }

    #[test]
    fn test_locus_name_extraction() {
        let gb = "LOCUS       MySeq                100 bp    DNA     linear   BCT\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("MySeq"));
    }

    #[test]
    fn test_circular_topology() {
        let gb = "LOCUS       CircSeq              100 bp    DNA     circular BCT\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert!(records[0].is_circular);
    }

    #[test]
    fn test_linear_topology() {
        let gb = "LOCUS       LinSeq               100 bp    DNA     linear   BCT\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert!(!records[0].is_circular);
    }

    #[test]
    fn test_definition_simple() {
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nDEFINITION  My cool sequence.\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].definition.as_deref(), Some("My cool sequence."));
    }

    #[test]
    fn test_definition_with_continuation() {
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nDEFINITION  A very long definition that\n            continues on the next line.\nACCESSION   X12345\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(
            records[0].definition.as_deref(),
            Some("A very long definition that continues on the next line.")
        );
        assert_eq!(records[0].accession.as_deref(), Some("X12345"));
    }

    #[test]
    fn test_accession_parsing() {
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nACCESSION   AB123456\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].accession.as_deref(), Some("AB123456"));
    }

    #[test]
    fn test_origin_sequence() {
        let gb = "LOCUS       S1 60 bp DNA linear UNK\nORIGIN\n        1 acgtacgtac gtacgtacgt acgtacgtac acgtacgtac gtacgtacgt acgtacgtac\n       61 ttttgggg\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].seq.len(), 68);
        assert!(records[0].seq.starts_with(b"acgtacgtac"));
        assert!(records[0].seq.ends_with(b"ttttgggg"));
    }

    #[test]
    fn test_multi_record() {
        let gb = "\
LOCUS       Seq1 4 bp DNA linear UNK
ORIGIN
        1 aaaa
//
LOCUS       Seq2 4 bp DNA circular UNK
ORIGIN
        1 cccc
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("Seq1"));
        assert_eq!(records[0].seq, b"aaaa");
        assert!(!records[0].is_circular);
        assert_eq!(records[1].name.as_deref(), Some("Seq2"));
        assert_eq!(records[1].seq, b"cccc");
        assert!(records[1].is_circular);
    }

    #[test]
    fn test_no_origin_empty_sequence() {
        let gb = "LOCUS       Empty 0 bp DNA linear UNK\nDEFINITION  No sequence.\n//\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert!(records[0].seq.is_empty());
    }

    #[test]
    fn test_minimal_fields() {
        let gb = "LOCUS       Min 4 bp DNA linear UNK\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("Min"));
        assert!(records[0].accession.is_none());
        assert!(records[0].definition.is_none());
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_features_skipped() {
        let gb = "\
LOCUS       Featured 8 bp DNA linear UNK
DEFINITION  Has features.
ACCESSION   F001
FEATURES             Location/Qualifiers
     gene            1..8
                     /locus_tag=\"test\"
     CDS             complement(1..8)
                     /product=\"hypothetical protein\"
                     /translation=\"MK\"
ORIGIN
        1 acgtacgt
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("Featured"));
        assert_eq!(records[0].definition.as_deref(), Some("Has features."));
        assert_eq!(records[0].accession.as_deref(), Some("F001"));
        assert_eq!(records[0].seq, b"acgtacgt");
    }

    #[test]
    fn test_empty_input() {
        let records = parse_all("");
        assert!(records.is_empty());
    }

    #[test]
    fn test_accession_multiple_tokens() {
        // ACCESSION line can have multiple accession numbers; we take the first
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nACCESSION   AB123456 AB789012\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].accession.as_deref(), Some("AB123456"));
    }

    #[test]
    fn test_realistic_full_record() {
        let gb = "\
LOCUS       AB012345             2400 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli strain K-12 16S ribosomal RNA gene, partial
            sequence.
ACCESSION   AB012345
VERSION     AB012345.1
DBLINK      BioProject: PRJNA12345
KEYWORDS    16S rRNA.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;
            Enterobacteriaceae; Escherichia.
REFERENCE   1  (bases 1 to 2400)
  AUTHORS   Smith,J. and Doe,A.
  TITLE     Some Paper Title
  JOURNAL   J. Bacteriol. 200 (1), 1-10 (2020)
   PUBMED   12345678
COMMENT     This is a comment that spans
            multiple lines.
FEATURES             Location/Qualifiers
     source          1..2400
                     /organism=\"Escherichia coli\"
                     /mol_type=\"genomic DNA\"
                     /strain=\"K-12\"
     gene            1..2400
                     /gene=\"16S rRNA\"
     rRNA            1..2400
                     /gene=\"16S rRNA\"
                     /product=\"16S ribosomal RNA\"
ORIGIN
        1 acgtacgtac gtacgtacgt
       21 ttttggggaa aacccctttt
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("AB012345"));
        assert_eq!(
            records[0].definition.as_deref(),
            Some("Escherichia coli strain K-12 16S ribosomal RNA gene, partial sequence.")
        );
        assert_eq!(records[0].accession.as_deref(), Some("AB012345"));
        assert!(records[0].is_circular);
        assert_eq!(records[0].seq.len(), 40);
        assert_eq!(
            std::str::from_utf8(&records[0].seq).unwrap(),
            "acgtacgtacgtacgtacgtttttggggaaaacccctttt"
        );
    }

    #[test]
    fn test_definition_followed_by_origin() {
        // DEFINITION with no ACCESSION between it and ORIGIN
        let gb = "\
LOCUS       S1 4 bp DNA linear UNK
DEFINITION  Direct to origin.
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(records[0].definition.as_deref(), Some("Direct to origin."));
        assert!(records[0].accession.is_none());
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_definition_followed_by_features() {
        // DEFINITION followed by FEATURES (not ACCESSION), ACCESSION appears later
        let gb = "\
LOCUS       S1 4 bp DNA linear UNK
DEFINITION  Has features first.
FEATURES             Location/Qualifiers
     gene            1..4
ACCESSION   LATE001
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(
            records[0].definition.as_deref(),
            Some("Has features first.")
        );
        // ACCESSION comes after FEATURES — should still be parsed
        assert_eq!(records[0].accession.as_deref(), Some("LATE001"));
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_definition_followed_by_record_terminator() {
        // DEFINITION continuation reads `//` as the next line
        let gb = "\
LOCUS       S1 0 bp DNA linear UNK
DEFINITION  No origin.
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].definition.as_deref(), Some("No origin."));
        assert!(records[0].seq.is_empty());
    }

    #[test]
    fn test_record_terminated_by_eof() {
        // No `//` at end — EOF terminates
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nORIGIN\n        1 acgt\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("S1"));
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_consecutive_locus_no_terminator() {
        // Second record starts without `//` terminating the first
        let gb = "\
LOCUS       First 4 bp DNA linear UNK
ORIGIN
        1 aaaa
LOCUS       Second 4 bp DNA linear UNK
ORIGIN
        1 cccc
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("First"));
        assert_eq!(records[0].seq, b"aaaa");
        assert_eq!(records[1].name.as_deref(), Some("Second"));
        assert_eq!(records[1].seq, b"cccc");
    }

    #[test]
    fn test_empty_accession_line() {
        // ACCESSION with no value
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nACCESSION   \nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert!(records[0].accession.is_none());
    }

    #[test]
    fn test_locus_no_topology_keyword() {
        // LOCUS line with neither "circular" nor "linear"
        let gb = "LOCUS       S1 4 bp DNA UNK\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert!(!records[0].is_circular);
    }

    #[test]
    fn test_from_bufreader_constructor() {
        let data = b"LOCUS       BR 4 bp DNA linear UNK\nORIGIN\n        1 acgt\n//\n";
        let buf_reader = BufReader::new(Cursor::new(data.as_slice()));
        let records: Vec<GenbankRecord> = GenbankReader::from_bufreader(buf_reader)
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("BR"));
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_preamble_before_locus() {
        // Garbage text before the first LOCUS line should be skipped
        let gb = "\
Some random text
Another line
LOCUS       S1 4 bp DNA linear UNK
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("S1"));
    }

    #[test]
    fn test_mixed_case_sequence_preserved() {
        // Parser should preserve the case from the file — consumer uppercases
        let gb = "LOCUS       S1 4 bp DNA linear UNK\nORIGIN\n        1 AcGt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].seq, b"AcGt");
    }

    #[test]
    fn test_origin_then_immediate_terminator() {
        // ORIGIN with no sequence lines before //
        let gb = "LOCUS       S1 0 bp DNA linear UNK\nORIGIN\n//\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert!(records[0].seq.is_empty());
    }

    #[test]
    fn test_accession_before_definition() {
        // Non-standard order: ACCESSION before DEFINITION
        let gb = "\
LOCUS       S1 4 bp DNA linear UNK
ACCESSION   ACC999
DEFINITION  After accession.
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(records[0].accession.as_deref(), Some("ACC999"));
        assert_eq!(records[0].definition.as_deref(), Some("After accession."));
    }

    #[test]
    fn test_short_locus_line() {
        // LOCUS line shorter than 12 characters — name extraction should handle gracefully
        let gb = "LOCUS\nORIGIN\n        1 acgt\n//\n";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert!(records[0].name.is_none());
        assert_eq!(records[0].seq, b"acgt");
    }

    #[test]
    fn test_three_records() {
        let gb = "\
LOCUS       A 2 bp DNA linear UNK
ORIGIN
        1 aa
//
LOCUS       B 2 bp DNA circular UNK
ORIGIN
        1 cc
//
LOCUS       C 2 bp DNA linear UNK
ORIGIN
        1 gg
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].name.as_deref(), Some("A"));
        assert_eq!(records[0].seq, b"aa");
        assert_eq!(records[1].name.as_deref(), Some("B"));
        assert_eq!(records[1].seq, b"cc");
        assert!(records[1].is_circular);
        assert_eq!(records[2].name.as_deref(), Some("C"));
        assert_eq!(records[2].seq, b"gg");
    }

    #[test]
    fn test_definition_three_continuation_lines() {
        let gb = "\
LOCUS       S1 4 bp DNA linear UNK
DEFINITION  Line one of the definition
            continues here on line two
            and finishes on line three.
ACCESSION   X999
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(
            records[0].definition.as_deref(),
            Some(
                "Line one of the definition continues here on line two and finishes on line three."
            )
        );
        assert_eq!(records[0].accession.as_deref(), Some("X999"));
    }

    #[test]
    fn test_large_sequence_many_origin_lines() {
        use std::fmt::Write;
        // Build a GenBank record with 600bp of sequence (10 ORIGIN lines)
        let seq_chunk = "acgtacgtac"; // 10 bases
        let mut gb = String::from("LOCUS       Big 600 bp DNA linear UNK\nORIGIN\n");
        let mut pos = 1;
        for _ in 0..10 {
            writeln!(
                gb,
                "{pos:>9} {seq_chunk} {seq_chunk} {seq_chunk} {seq_chunk} {seq_chunk} {seq_chunk}"
            )
            .unwrap();
            pos += 60;
        }
        gb.push_str("//\n");

        let records = parse_all(&gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq.len(), 600);
        // All bases should be alphabetic
        assert!(records[0].seq.iter().all(u8::is_ascii_alphabetic));
    }

    #[test]
    fn test_sequence_with_n_characters() {
        // GenBank sequences can contain N for unknown bases
        let gb = "LOCUS       S1 8 bp DNA linear UNK\nORIGIN\n        1 acgtnnnn\n//\n";
        let records = parse_all(gb);
        assert_eq!(records[0].seq, b"acgtnnnn");
    }

    #[test]
    fn test_definition_eof_during_continuation() {
        // EOF during DEFINITION continuation (no terminator at all)
        let gb =
            "LOCUS       S1 0 bp DNA linear UNK\nDEFINITION  Truncated at\n            end of file";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].definition.as_deref(),
            Some("Truncated at end of file")
        );
        assert!(records[0].seq.is_empty());
    }

    #[test]
    fn test_only_whitespace_between_records() {
        let gb = "\
LOCUS       R1 2 bp DNA linear UNK
ORIGIN
        1 aa
//

LOCUS       R2 2 bp DNA linear UNK
ORIGIN
        1 cc
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("R1"));
        assert_eq!(records[1].name.as_deref(), Some("R2"));
    }

    #[test]
    fn test_max_line_length_enforced() {
        // A GenBank "file" with a line exceeding MAX_LINE_LENGTH
        let long_line = "A".repeat(MAX_LINE_LENGTH + 1);
        let gb = format!("LOCUS       S1 4 bp DNA linear UNK\n{long_line}\n//\n");
        let reader = GenbankReader::new(Cursor::new(gb.as_bytes()));
        let results: Vec<_> = reader.collect();
        // Should have an error for the long line
        assert!(
            results.iter().any(std::result::Result::is_err),
            "Should reject line exceeding MAX_LINE_LENGTH"
        );
    }

    #[test]
    fn test_max_line_length_ok_at_limit() {
        // Line exactly at MAX_LINE_LENGTH should be OK
        let line = "A".repeat(MAX_LINE_LENGTH);
        let gb = format!("LOCUS       S1 4 bp DNA linear UNK\n{line}\n//\n");
        let reader = GenbankReader::new(Cursor::new(gb.as_bytes()));
        let results: Vec<_> = reader.collect();
        // All should be Ok (no error)
        assert!(
            results.iter().all(std::result::Result::is_ok),
            "Line at exactly MAX_LINE_LENGTH should be accepted"
        );
    }

    // ========================================================================
    // Slice-based fast path tests
    // ========================================================================

    #[test]
    fn test_slice_basic() {
        let gb = b"LOCUS       MySeq                100 bp    DNA     linear   BCT\nDEFINITION  A test.\nACCESSION   X123\nORIGIN\n        1 acgtacgt\n//\n";
        let records = parse_genbank_slice(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("MySeq"));
        assert_eq!(records[0].definition.as_deref(), Some("A test."));
        assert_eq!(records[0].accession.as_deref(), Some("X123"));
        assert!(!records[0].is_circular);
        assert_eq!(records[0].seq, b"acgtacgt");
    }

    #[test]
    fn test_slice_multi_record() {
        let gb = b"\
LOCUS       Seq1 4 bp DNA linear UNK
ORIGIN
        1 aaaa
//
LOCUS       Seq2 4 bp DNA circular UNK
ORIGIN
        1 cccc
//
";
        let records = parse_genbank_slice(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("Seq1"));
        assert_eq!(records[0].seq, b"aaaa");
        assert!(!records[0].is_circular);
        assert_eq!(records[1].name.as_deref(), Some("Seq2"));
        assert_eq!(records[1].seq, b"cccc");
        assert!(records[1].is_circular);
    }

    #[test]
    fn test_slice_features_skipped() {
        let gb = b"\
LOCUS       Featured 8 bp DNA linear UNK
DEFINITION  Has features.
ACCESSION   F001
FEATURES             Location/Qualifiers
     gene            1..8
                     /locus_tag=\"test\"
     CDS             complement(1..8)
                     /product=\"hypothetical protein\"
                     /translation=\"MK\"
ORIGIN
        1 acgtacgt
//
";
        let records = parse_genbank_slice(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("Featured"));
        assert_eq!(records[0].definition.as_deref(), Some("Has features."));
        assert_eq!(records[0].accession.as_deref(), Some("F001"));
        assert_eq!(records[0].seq, b"acgtacgt");
    }

    #[test]
    fn test_slice_consecutive_no_terminator() {
        let gb = b"\
LOCUS       First 4 bp DNA linear UNK
ORIGIN
        1 aaaa
LOCUS       Second 4 bp DNA linear UNK
ORIGIN
        1 cccc
//
";
        let records = parse_genbank_slice(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("First"));
        assert_eq!(records[0].seq, b"aaaa");
        assert_eq!(records[1].name.as_deref(), Some("Second"));
        assert_eq!(records[1].seq, b"cccc");
    }

    #[test]
    fn test_slice_parity_with_streaming() {
        let gb = "\
LOCUS       AB012345             2400 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli strain K-12 16S ribosomal RNA gene, partial
            sequence.
ACCESSION   AB012345
VERSION     AB012345.1
FEATURES             Location/Qualifiers
     source          1..2400
                     /organism=\"Escherichia coli\"
     gene            1..2400
                     /gene=\"16S rRNA\"
ORIGIN
        1 acgtacgtac gtacgtacgt
       21 ttttggggaa aacccctttt
//
";
        // Streaming path
        let streaming: Vec<GenbankRecord> = GenbankReader::new(Cursor::new(gb.as_bytes()))
            .collect::<Result<Vec<_>>>()
            .unwrap();
        // Slice path
        let sliced = parse_genbank_slice(gb.as_bytes());

        assert_eq!(streaming.len(), sliced.len());
        for (s, sl) in streaming.iter().zip(sliced.iter()) {
            assert_eq!(s.name, sl.name);
            assert_eq!(s.accession, sl.accession);
            assert_eq!(s.definition, sl.definition);
            assert_eq!(s.is_circular, sl.is_circular);
            assert_eq!(s.seq, sl.seq);
        }
    }

    #[test]
    fn test_slice_realistic_full_record() {
        let gb = b"\
LOCUS       AB012345             2400 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli strain K-12 16S ribosomal RNA gene, partial
            sequence.
ACCESSION   AB012345
VERSION     AB012345.1
DBLINK      BioProject: PRJNA12345
KEYWORDS    16S rRNA.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;
            Enterobacteriaceae; Escherichia.
REFERENCE   1  (bases 1 to 2400)
  AUTHORS   Smith,J. and Doe,A.
  TITLE     Some Paper Title
  JOURNAL   J. Bacteriol. 200 (1), 1-10 (2020)
   PUBMED   12345678
COMMENT     This is a comment that spans
            multiple lines.
FEATURES             Location/Qualifiers
     source          1..2400
                     /organism=\"Escherichia coli\"
                     /mol_type=\"genomic DNA\"
                     /strain=\"K-12\"
     gene            1..2400
                     /gene=\"16S rRNA\"
     rRNA            1..2400
                     /gene=\"16S rRNA\"
                     /product=\"16S ribosomal RNA\"
ORIGIN
        1 acgtacgtac gtacgtacgt
       21 ttttggggaa aacccctttt
//
";
        let records = parse_genbank_slice(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name.as_deref(), Some("AB012345"));
        assert_eq!(
            records[0].definition.as_deref(),
            Some("Escherichia coli strain K-12 16S ribosomal RNA gene, partial sequence.")
        );
        assert_eq!(records[0].accession.as_deref(), Some("AB012345"));
        assert!(records[0].is_circular);
        assert_eq!(records[0].seq.len(), 40);
        assert_eq!(
            std::str::from_utf8(&records[0].seq).unwrap(),
            "acgtacgtacgtacgtacgtttttggggaaaacccctttt"
        );
    }

    // ========================================================================
    // Parity tests: verify streaming and slice parsers produce identical output
    // ========================================================================

    /// Assert that streaming and slice parsers produce field-by-field identical
    /// records for the given input bytes.
    fn assert_parser_parity(input: &[u8]) {
        // Streaming path
        let streaming: Vec<GenbankRecord> = GenbankReader::new(Cursor::new(input))
            .collect::<Result<Vec<_>>>()
            .expect("streaming parser should not error");
        // Slice path
        let sliced = parse_genbank_slice(input);

        assert_eq!(
            streaming.len(),
            sliced.len(),
            "Record count mismatch: streaming={}, slice={}",
            streaming.len(),
            sliced.len()
        );

        for (i, (s, sl)) in streaming.iter().zip(sliced.iter()).enumerate() {
            assert_eq!(s.name, sl.name, "Record {i}: name mismatch");
            assert_eq!(s.accession, sl.accession, "Record {i}: accession mismatch");
            assert_eq!(
                s.definition, sl.definition,
                "Record {i}: definition mismatch"
            );
            assert_eq!(
                s.is_circular, sl.is_circular,
                "Record {i}: circular mismatch"
            );
            assert_eq!(s.seq, sl.seq, "Record {i}: sequence mismatch");
        }
    }

    #[test]
    fn parity_realistic_record() {
        let input = b"\
LOCUS       AB012345             2400 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  Escherichia coli strain K-12 16S ribosomal RNA gene, partial
            sequence.
ACCESSION   AB012345
VERSION     AB012345.1
DBLINK      BioProject: PRJNA12345
KEYWORDS    16S rRNA.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;
            Enterobacteriaceae; Escherichia.
REFERENCE   1  (bases 1 to 2400)
  AUTHORS   Smith,J. and Doe,A.
  TITLE     Some Paper Title
  JOURNAL   J. Bacteriol. 200 (1), 1-10 (2020)
   PUBMED   12345678
COMMENT     This is a comment that spans
            multiple lines.
FEATURES             Location/Qualifiers
     source          1..2400
                     /organism=\"Escherichia coli\"
                     /mol_type=\"genomic DNA\"
                     /strain=\"K-12\"
     gene            1..2400
                     /gene=\"16S rRNA\"
     rRNA            1..2400
                     /gene=\"16S rRNA\"
                     /product=\"16S ribosomal RNA\"
ORIGIN
        1 acgtacgtac gtacgtacgt
       21 ttttggggaa aacccctttt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_multi_record() {
        let input = b"\
LOCUS       Linear1 4 bp DNA linear UNK
DEFINITION  First linear record.
ACCESSION   LIN001
ORIGIN
        1 aaaa
//
LOCUS       Circ1 4 bp DNA circular UNK
DEFINITION  Second circular record.
ACCESSION   CIR001
ORIGIN
        1 cccc
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_no_origin() {
        let input = b"\
LOCUS       NoOrig 0 bp DNA linear UNK
DEFINITION  Record without ORIGIN section.
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_truncated_eof() {
        // Record terminated by EOF instead of //
        let input = b"\
LOCUS       TruncEOF 4 bp DNA linear UNK
DEFINITION  Truncated record.
ACCESSION   TRUNC001
ORIGIN
        1 acgt
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_missing_terminator() {
        // Two consecutive LOCUS lines without // between them
        let input = b"\
LOCUS       First 4 bp DNA linear UNK
ORIGIN
        1 aaaa
LOCUS       Second 4 bp DNA circular UNK
ORIGIN
        1 cccc
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_definition_continuation() {
        // DEFINITION spanning 3+ continuation lines
        let input = b"\
LOCUS       DefCont 4 bp DNA linear UNK
DEFINITION  This is a very long definition that spans
            multiple continuation lines in the GenBank
            format, which is quite common for detailed
            sequence descriptions.
ACCESSION   DEF001
ORIGIN
        1 acgt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_empty_input() {
        assert_parser_parity(b"");
    }

    #[test]
    fn parity_minimal_record() {
        // Just LOCUS + ORIGIN + //
        let input = b"\
LOCUS       MinRec 4 bp DNA linear UNK
ORIGIN
        1 acgt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_features_between_header_and_origin() {
        // FEATURES section with gene/CDS annotations between ACCESSION and ORIGIN
        let input = b"\
LOCUS       FeatRec 8 bp DNA linear UNK
DEFINITION  Record with features.
ACCESSION   FEAT001
FEATURES             Location/Qualifiers
     gene            1..8
                     /gene=\"testGene\"
     CDS             1..6
                     /product=\"hypothetical protein\"
                     /translation=\"MK\"
ORIGIN
        1 acgtacgt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn parity_accession_before_definition() {
        // Non-standard field order
        let input = b"\
LOCUS       RevOrder 4 bp DNA linear UNK
ACCESSION   REV001
DEFINITION  Definition after accession.
ORIGIN
        1 acgt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn test_definition_then_locus_no_terminator() {
        // DEFINITION continuation stops at LOCUS (no ORIGIN, no //)
        // This exercises the `has_pending_locus` branch after DEFINITION parsing
        let gb = "\
LOCUS       First 0 bp DNA linear UNK
DEFINITION  First record definition.
LOCUS       Second 4 bp DNA linear UNK
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name.as_deref(), Some("First"));
        assert_eq!(
            records[0].definition.as_deref(),
            Some("First record definition.")
        );
        assert!(records[0].seq.is_empty());
        assert_eq!(records[1].name.as_deref(), Some("Second"));
        assert_eq!(records[1].seq, b"acgt");
    }

    #[test]
    fn parity_definition_then_locus_no_terminator() {
        let input = b"\
LOCUS       First 0 bp DNA linear UNK
DEFINITION  First record definition.
LOCUS       Second 4 bp DNA linear UNK
ORIGIN
        1 acgt
//
";
        assert_parser_parity(input);
    }

    #[test]
    fn test_definition_multiline_then_origin() {
        // Multi-line DEFINITION followed directly by ORIGIN (no ACCESSION)
        let gb = "\
LOCUS       S1 4 bp DNA linear UNK
DEFINITION  First line of definition
            second line of definition
            third line of definition.
ORIGIN
        1 acgt
//
";
        let records = parse_all(gb);
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].definition.as_deref(),
            Some("First line of definition second line of definition third line of definition.")
        );
        assert_eq!(records[0].seq, b"acgt");
    }
}
