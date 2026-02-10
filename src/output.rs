//! Output formatting for FASTA sequences and TSV statistics tables.
//!
//! Provides streaming FASTA writers (plain and gzip-compressed) and
//! tab-separated summary output with per-product alignment details.

use anyhow::{Context, Result};
use gzp::{ZBuilder, deflate::Bgzf};
use sassy::Strand;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::errors::ValidationError;
use crate::pcr::PcrProduct;

/// Line width for FASTA sequence wrapping
const FASTA_LINE_WIDTH: usize = 80;

/// Streaming FASTA writer (plain or gzipped)
pub enum FastaWriter {
    /// Uncompressed output through a buffered file writer.
    Plain(BufWriter<File>),
    /// Multi-threaded BGZF-compressed output via gzp.
    Gzip(Box<dyn gzp::ZWriter<File>>),
}

impl FastaWriter {
    /// Create a new FASTA writer for the given output path.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    pub fn new(output_path: &Path, threads: usize) -> Result<Self> {
        let is_gzipped = output_path
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

        let file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

        if is_gzipped {
            // Use BGZF format for parallel compression and future parallel decompression
            // BGZF is gzip-compatible (any gzip reader can decompress it)
            let writer: Box<dyn gzp::ZWriter<File>> = ZBuilder::<Bgzf, _>::new()
                .num_threads(threads)
                .from_writer(file);
            Ok(FastaWriter::Gzip(writer))
        } else {
            Ok(FastaWriter::Plain(BufWriter::with_capacity(
                64 * 1024,
                file,
            )))
        }
    }

    /// Write a single PCR product as a FASTA record.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the underlying stream fails.
    pub fn write_product(&mut self, product: &PcrProduct) -> Result<()> {
        match self {
            FastaWriter::Plain(writer) => write_fasta_record(writer, product),
            FastaWriter::Gzip(writer) => write_fasta_record(writer, product),
        }
    }

    /// Flush and finalize the writer.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or finalizing the output stream fails.
    pub fn finish(self) -> Result<()> {
        match self {
            FastaWriter::Plain(mut writer) => {
                writer.flush()?;
                Ok(())
            }
            FastaWriter::Gzip(mut writer) => {
                writer
                    .finish()
                    .with_context(|| "Failed to finalize gzip output")?;
                Ok(())
            }
        }
    }
}

/// Write PCR products to a FASTA file (with optional gzip compression).
///
/// # Errors
///
/// Returns an error if the output file cannot be created or written to.
pub fn write_fasta(products: &[PcrProduct], output_path: &Path, threads: usize) -> Result<()> {
    let mut writer = FastaWriter::new(output_path, threads)?;
    for product in products {
        writer.write_product(product)?;
    }
    writer.finish()?;
    log::info!(
        "Wrote {} products to {}",
        products.len(),
        output_path.display()
    );
    Ok(())
}

/// Write a single FASTA record.
///
/// # Errors
///
/// Returns an error if writing to the stream fails.
pub fn write_fasta_record<W: Write>(writer: &mut W, product: &PcrProduct) -> Result<()> {
    // Write header
    let strand_str = match product.strand {
        Strand::Fwd => "+",
        Strand::Rc => "-",
    };
    writeln!(
        writer,
        ">{}\tpos={}-{}\tstrand={}\tlen={}",
        product.header(),
        product.original_start,
        product.original_end,
        strand_str,
        product.full_length
    )?;

    // Write sequence with line wrapping
    for chunk in product.sequence.chunks(FASTA_LINE_WIDTH) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }

    Ok(())
}

/// Write PCR products to stdout as FASTA.
///
/// # Errors
///
/// Returns an error if writing to stdout fails.
pub fn write_fasta_stdout(products: &[PcrProduct]) -> Result<()> {
    let stdout = std::io::stdout();
    let mut writer = BufWriter::new(stdout.lock());

    for product in products {
        write_fasta_record(&mut writer, product)?;
    }

    writer.flush()?;
    Ok(())
}

/// Validate that an output path is writable before starting expensive processing.
///
/// Checks:
/// 1. Parent directory exists (returns `ValidationError::OutputDirMissing` if not).
/// 2. File can be created at the path (returns `ValidationError::OutputNotWritable` on failure).
///    A probe file is created and immediately removed.
///
/// # Errors
///
/// Returns a `ValidationError` if the output path is not writable.
pub fn validate_output_writable(path: &Path) -> Result<()> {
    // Check parent directory exists
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
        && !parent.exists()
    {
        return Err(ValidationError::OutputDirMissing {
            path: parent.to_path_buf(),
        }
        .into());
    }

    // Probe: try creating the file, then immediately remove it
    match File::create(path) {
        Ok(_) => {
            // Clean up the probe file
            let _ = std::fs::remove_file(path);
            Ok(())
        }
        Err(e) => Err(ValidationError::OutputNotWritable {
            path: path.to_path_buf(),
            source: e,
        }
        .into()),
    }
}

/// TSV output columns
const TSV_HEADER: &str = "amplicon_id\treference_id\tsource_file\tprimer_name\tproduct_len\tfull_len\t\
fwd_start\tfwd_end\tfwd_mismatches\tfwd_identity\tfwd_cigar\t\
rev_start\trev_end\trev_mismatches\trev_identity\trev_cigar\t\
strand\tis_circular_wrap\tproduct_seq";

/// Streaming TSV writer
pub struct TsvWriter {
    writer: BufWriter<File>,
}

impl TsvWriter {
    /// Create a new TSV writer for the given output path.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    pub fn new(output_path: &Path) -> Result<Self> {
        let file = File::create(output_path)
            .with_context(|| format!("Failed to create TSV file: {}", output_path.display()))?;
        let mut writer = BufWriter::with_capacity(64 * 1024, file);
        writeln!(writer, "{TSV_HEADER}")?;
        Ok(Self { writer })
    }

    /// Write a single PCR product as a TSV record.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the underlying stream fails.
    pub fn write_product(&mut self, product: &PcrProduct) -> Result<()> {
        write_tsv_record(&mut self.writer, product)
    }

    /// Flush and finalize the writer.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the output stream fails.
    pub fn finish(mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Write TSV statistics file.
///
/// # Errors
///
/// Returns an error if the output file cannot be created or written to.
pub fn write_tsv(products: &[PcrProduct], output_path: &Path) -> Result<()> {
    let mut writer = TsvWriter::new(output_path)?;
    for product in products {
        writer.write_product(product)?;
    }
    writer.finish()?;
    log::info!("Wrote TSV statistics to {}", output_path.display());
    Ok(())
}

/// Write a single TSV record
fn write_tsv_record<W: Write>(writer: &mut W, product: &PcrProduct) -> Result<()> {
    let strand_str = match product.strand {
        Strand::Fwd => "+",
        Strand::Rc => "-",
    };

    let seq_str = String::from_utf8_lossy(&product.sequence);

    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
        product.header(),
        product.reference_header,
        product.source_file,
        product.primer_name,
        product.len(),
        product.full_length,
        product.fwd_match.start,
        product.fwd_match.end,
        product.fwd_match.edit_distance,
        product.fwd_match.identity,
        product.fwd_match.cigar,
        product.rev_match.start,
        product.rev_match.end,
        product.rev_match.edit_distance,
        product.rev_match.identity,
        product.rev_match.cigar,
        strand_str,
        product.is_circular_wrap,
        seq_str,
    )?;

    Ok(())
}

/// Summary statistics for a completed PCR analysis run.
#[derive(Debug, Default)]
pub struct RunSummary {
    /// Number of input sequences (FASTA records or `GenBank` entries) processed.
    pub total_sequences: usize,
    /// Number of primer pairs used in the search.
    pub total_primers: usize,
    /// Total PCR products found across all sequences and primers.
    pub total_products: usize,
    /// Product count grouped by primer pair name.
    pub products_per_primer: Vec<(String, usize)>,
    /// Product count grouped by reference sequence header.
    pub products_per_reference: Vec<(String, usize)>,
}

impl RunSummary {
    /// Create summary from pre-aggregated counts
    #[must_use]
    pub fn from_counts(
        total_sequences: usize,
        total_primers: usize,
        total_products: usize,
        primer_counts: std::collections::HashMap<String, usize>,
        ref_counts: std::collections::HashMap<String, usize>,
    ) -> Self {
        let mut products_per_primer: Vec<_> = primer_counts.into_iter().collect();
        products_per_primer.sort_by(|a, b| b.1.cmp(&a.1));

        let mut products_per_reference: Vec<_> = ref_counts.into_iter().collect();
        products_per_reference.sort_by(|a, b| b.1.cmp(&a.1));

        Self {
            total_sequences,
            total_primers,
            total_products,
            products_per_primer,
            products_per_reference,
        }
    }

    /// Create summary from products
    #[must_use]
    pub fn from_products(
        products: &[PcrProduct],
        num_sequences: usize,
        num_primers: usize,
    ) -> Self {
        use std::collections::HashMap;

        let mut primer_counts: HashMap<String, usize> = HashMap::new();
        let mut ref_counts: HashMap<String, usize> = HashMap::new();

        for product in products {
            *primer_counts
                .entry(product.primer_name.clone())
                .or_default() += 1;
            *ref_counts
                .entry(product.reference_header.clone())
                .or_default() += 1;
        }

        let mut products_per_primer: Vec<_> = primer_counts.into_iter().collect();
        products_per_primer.sort_by(|a, b| b.1.cmp(&a.1));

        let mut products_per_reference: Vec<_> = ref_counts.into_iter().collect();
        products_per_reference.sort_by(|a, b| b.1.cmp(&a.1));

        Self {
            total_sequences: num_sequences,
            total_primers: num_primers,
            total_products: products.len(),
            products_per_primer,
            products_per_reference,
        }
    }

    /// Print summary to log (requires -v flag)
    pub fn log_summary(&self) {
        log::info!("=== Run Summary ===");
        log::info!("Input sequences: {}", self.total_sequences);
        log::info!("Primer pairs: {}", self.total_primers);
        log::info!("PCR products found: {}", self.total_products);

        if !self.products_per_primer.is_empty() {
            log::info!("Products by primer:");
            for (primer, count) in &self.products_per_primer {
                log::info!("  {primer}: {count}");
            }
        }

        if self.products_per_reference.len() <= 10 {
            log::info!("Products by reference:");
            for (ref_name, count) in &self.products_per_reference {
                log::info!("  {ref_name}: {count}");
            }
        } else {
            log::info!(
                "Products found in {} different references",
                self.products_per_reference.len()
            );
        }
    }

    /// Print summary to stderr (always, regardless of log level)
    pub fn print_stderr(&self) {
        eprintln!();
        eprintln!("=== Amplirust Summary ===");
        eprintln!("Input sequences:  {}", self.total_sequences);
        eprintln!("Primer pairs:     {}", self.total_primers);
        eprintln!("Products found:   {}", self.total_products);

        if !self.products_per_primer.is_empty() {
            eprintln!();
            eprintln!("Products by primer:");
            for (primer, count) in &self.products_per_primer {
                eprintln!("  {primer}: {count}");
            }
        }

        if !self.products_per_reference.is_empty() {
            eprintln!();
            if self.products_per_reference.len() <= 10 {
                eprintln!("Products by reference:");
                for (ref_name, count) in &self.products_per_reference {
                    eprintln!("  {ref_name}: {count}");
                }
            } else {
                eprintln!(
                    "Products found in {} different references",
                    self.products_per_reference.len()
                );
            }
        }
        eprintln!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matcher::PrimerMatch;
    use tempfile::NamedTempFile;

    fn make_test_product(case: usize) -> PcrProduct {
        PcrProduct {
            reference_header: "test_ref".to_string(),
            source_file: "test.fa".to_string(),
            primer_name: "test_primer".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
            full_length: 12,
            fwd_match: PrimerMatch {
                start: 0,
                end: 4,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            rev_match: PrimerMatch {
                start: 8,
                end: 12,
                edit_distance: 0,
                strand: Strand::Fwd,
                cigar: "4=".to_string(),
                identity: 1.0,
            },
            strand: Strand::Fwd,
            is_circular_wrap: false,
            original_start: 0,
            original_end: 12,
            case_number: case,
        }
    }

    #[test]
    fn test_write_fasta_plain() {
        let products = vec![make_test_product(1), make_test_product(2)];
        let temp = NamedTempFile::with_suffix(".fasta").unwrap();

        write_fasta(&products, temp.path(), 1).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains(">test_ref:test_primer:1"));
        assert!(content.contains(">test_ref:test_primer:2"));
        assert!(content.contains("ACGTACGTACGT"));
    }

    #[test]
    fn test_write_tsv() {
        let products = vec![make_test_product(1)];
        let temp = NamedTempFile::with_suffix(".tsv").unwrap();

        write_tsv(&products, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("amplicon_id"));
        assert!(content.contains("test_ref:test_primer:1"));
    }

    #[test]
    fn test_run_summary() {
        let products = vec![make_test_product(1), make_test_product(2)];

        let summary = RunSummary::from_products(&products, 10, 2);
        assert_eq!(summary.total_products, 2);
        assert_eq!(summary.total_sequences, 10);
    }

    #[test]
    fn test_tsv_all_columns_present() {
        // Verify all 19 columns in header and data row
        let product = make_test_product(1);
        let mut output = Vec::new();

        // Write header manually for test
        writeln!(output, "{TSV_HEADER}").unwrap();
        write_tsv_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines.len(), 2, "Should have header and one data row");

        // Check header has 19 columns
        let header_cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(
            header_cols.len(),
            19,
            "Header should have 19 columns, got {}",
            header_cols.len()
        );

        // Check data row has 19 columns
        let data_cols: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(
            data_cols.len(),
            19,
            "Data row should have 19 columns, got {}",
            data_cols.len()
        );

        // Verify specific column names are present
        assert!(header_cols.contains(&"amplicon_id"));
        assert!(header_cols.contains(&"reference_id"));
        assert!(header_cols.contains(&"source_file"));
        assert!(header_cols.contains(&"primer_name"));
        assert!(header_cols.contains(&"product_len"));
        assert!(header_cols.contains(&"full_len"));
        assert!(header_cols.contains(&"fwd_start"));
        assert!(header_cols.contains(&"fwd_end"));
        assert!(header_cols.contains(&"fwd_mismatches"));
        assert!(header_cols.contains(&"fwd_identity"));
        assert!(header_cols.contains(&"fwd_cigar"));
        assert!(header_cols.contains(&"rev_start"));
        assert!(header_cols.contains(&"rev_end"));
        assert!(header_cols.contains(&"rev_mismatches"));
        assert!(header_cols.contains(&"rev_identity"));
        assert!(header_cols.contains(&"rev_cigar"));
        assert!(header_cols.contains(&"strand"));
        assert!(header_cols.contains(&"is_circular_wrap"));
        assert!(header_cols.contains(&"product_seq"));
    }

    #[test]
    fn test_tsv_special_characters_in_header() {
        // Test reference header with quotes and spaces (tabs would break TSV format)
        let mut product = make_test_product(1);
        product.reference_header = "ref with spaces and \"quotes\" here".to_string();

        let mut output = Vec::new();
        write_tsv_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();

        // The TSV should still be parseable (tabs separate columns)
        let cols: Vec<&str> = content.trim().split('\t').collect();
        assert_eq!(cols.len(), 19, "Should still have 19 columns");

        // The reference_id column should contain the header
        assert!(content.contains("ref with spaces"));
        assert!(content.contains("quotes"));
    }

    #[test]
    fn test_tsv_tab_in_header_breaks_format() {
        // Document that tabs in headers break TSV format (expected behavior)
        let mut product = make_test_product(1);
        product.reference_header = "ref with\ttab".to_string();

        let mut output = Vec::new();
        write_tsv_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();

        // Tab in header creates extra column - this is known/expected behavior
        let cols: Vec<&str> = content.trim().split('\t').collect();
        assert!(
            cols.len() > 19,
            "Tab in header breaks TSV format (creates extra columns)"
        );
    }

    #[test]
    fn test_validate_output_nonexistent_dir() {
        let result = validate_output_writable(Path::new("/nonexistent/dir/output.fasta"));
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("does not exist"),
            "Expected 'does not exist' in error: {err}"
        );
    }

    #[test]
    fn test_validate_output_writable_ok() {
        let temp = tempfile::TempDir::new().unwrap();
        let path = temp.path().join("test_output.fasta");
        let result = validate_output_writable(&path);
        assert!(result.is_ok());
        // Probe file should have been cleaned up
        assert!(!path.exists());
    }

    #[test]
    fn test_fasta_line_wrapping() {
        // Sequence > 80bp should wrap correctly
        let mut product = make_test_product(1);
        // Create a sequence longer than 80bp
        product.sequence = b"ACGT".repeat(25); // 100bp

        let mut output = Vec::new();
        write_fasta_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        // Should have header + at least 2 sequence lines (100bp / 80 = 2 lines)
        assert!(
            lines.len() >= 3,
            "Should have header + wrapped sequence lines"
        );

        // First line is header
        assert!(lines[0].starts_with('>'));

        // Sequence lines should be max 80 chars
        for line in &lines[1..] {
            assert!(
                line.len() <= 80,
                "Sequence line should be max 80 chars, got {}",
                line.len()
            );
        }

        // Total sequence length should match
        let seq_content: String = lines[1..].join("");
        assert_eq!(seq_content.len(), 100);
    }

    #[test]
    fn test_fasta_short_sequence_no_wrap() {
        // Sequence < 80bp should be on single line
        let product = make_test_product(1); // 12bp sequence

        let mut output = Vec::new();
        write_fasta_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        // Should have header + 1 sequence line
        assert_eq!(lines.len(), 2, "Short sequence should be on one line");
        assert!(lines[0].starts_with('>'));
        assert_eq!(lines[1], "ACGTACGTACGT");
    }

    #[test]
    fn test_fasta_header_format() {
        let product = make_test_product(1);

        let mut output = Vec::new();
        write_fasta_record(&mut output, &product).unwrap();

        let content = String::from_utf8(output).unwrap();
        let header_line = content.lines().next().unwrap();

        // Header should contain expected fields
        assert!(header_line.starts_with('>'));
        assert!(header_line.contains("pos="));
        assert!(header_line.contains("strand="));
        assert!(header_line.contains("len="));
    }
}
