use anyhow::{Context, Result};
use gzp::{ZBuilder, deflate::Gzip};
use sassy::Strand;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::pcr::PcrProduct;

/// Line width for FASTA sequence wrapping
const FASTA_LINE_WIDTH: usize = 80;

/// Streaming FASTA writer (plain or gzipped)
pub enum FastaWriter {
    Plain(BufWriter<File>),
    Gzip(Box<dyn gzp::ZWriter>),
}

impl FastaWriter {
    pub fn new(output_path: &Path, threads: usize) -> Result<Self> {
        let is_gzipped = output_path
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

        let file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

        if is_gzipped {
            let writer: Box<dyn gzp::ZWriter> = ZBuilder::<Gzip, _>::new()
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

    pub fn write_product(&mut self, product: &PcrProduct) -> Result<()> {
        match self {
            FastaWriter::Plain(writer) => write_fasta_record(writer, product),
            FastaWriter::Gzip(writer) => write_fasta_record(writer, product),
        }
    }

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

/// Write PCR products to a FASTA file (with optional gzip compression)
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

/// Write a single FASTA record
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

/// Write PCR products to stdout as FASTA
pub fn write_fasta_stdout(products: &[PcrProduct]) -> Result<()> {
    let stdout = std::io::stdout();
    let mut writer = BufWriter::new(stdout.lock());

    for product in products {
        write_fasta_record(&mut writer, product)?;
    }

    writer.flush()?;
    Ok(())
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
    pub fn new(output_path: &Path) -> Result<Self> {
        let file = File::create(output_path)
            .with_context(|| format!("Failed to create TSV file: {}", output_path.display()))?;
        let mut writer = BufWriter::with_capacity(64 * 1024, file);
        writeln!(writer, "{}", TSV_HEADER)?;
        Ok(Self { writer })
    }

    pub fn write_product(&mut self, product: &PcrProduct) -> Result<()> {
        write_tsv_record(&mut self.writer, product)
    }

    pub fn finish(mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Write TSV statistics file
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

    let seq_str = std::str::from_utf8(&product.sequence).unwrap_or("<invalid UTF-8>");

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

/// Summary statistics for a run
#[derive(Debug, Default)]
pub struct RunSummary {
    pub total_sequences: usize,
    pub total_primers: usize,
    pub total_products: usize,
    pub products_per_primer: Vec<(String, usize)>,
    pub products_per_reference: Vec<(String, usize)>,
}

impl RunSummary {
    /// Create summary from pre-aggregated counts
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
                log::info!("  {}: {}", primer, count);
            }
        }

        if self.products_per_reference.len() <= 10 {
            log::info!("Products by reference:");
            for (ref_name, count) in &self.products_per_reference {
                log::info!("  {}: {}", ref_name, count);
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
                eprintln!("  {}: {}", primer, count);
            }
        }

        if !self.products_per_reference.is_empty() {
            eprintln!();
            if self.products_per_reference.len() <= 10 {
                eprintln!("Products by reference:");
                for (ref_name, count) in &self.products_per_reference {
                    eprintln!("  {}: {}", ref_name, count);
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
}
