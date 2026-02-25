//! Command-line argument parsing for amplirust.
//!
//! Uses clap derive macros for argument definition and validation.

use clap::Parser;
use std::path::PathBuf;

/// Amplirust: In-silico PCR primer matching and product extraction tool
#[derive(Parser, Debug, Clone)]
#[command(name = "amplirust")]
#[command(author, version, about, long_about = None)]
pub struct Args {
    // ==================== INPUT OPTIONS ====================
    /// Input sequence files (comma-separated, glob patterns supported)
    /// Supported formats: FASTA (.fasta, .fa, .fna, .ffn, .fas) and
    /// `GenBank` (.gb, .gbk, .gbff, .genbank, .gbf). Gzip (.gz) supported.
    /// Unrecognized file types in glob results are skipped with a warning.
    #[arg(short, long, required = true, value_delimiter = ',')]
    pub input: Vec<String>,

    /// Primers as "name:forward:reverse" or path to CSV file
    /// CSV format: name,forward,reverse (with header)
    /// Example: "16S:AGAGTTTGATCMTGGCTCAG:TACGGYTACCTTGTTACGACTT"
    #[arg(short, long, required = true)]
    pub primers: String,

    /// Treat sequences as circular genomes (allows wrap-around products)
    #[arg(long, default_value_t = false)]
    pub circular: bool,

    /// Maximum decompressed file size in bytes (0 = unlimited).
    /// Prevents memory exhaustion from malicious compressed files.
    #[arg(long, default_value_t = 4_294_967_296)]
    pub max_decompression_size: u64,

    // ==================== MATCHING OPTIONS ====================
    /// Maximum edit distance (mismatches + indels) for primer matching
    #[arg(short = 'k', long, default_value_t = 2)]
    pub max_errors: usize,

    /// Minimum percentage identity for a primer match (0.0 - 1.0, 0 = disabled)
    #[arg(long, default_value_t = 0.0)]
    pub min_identity: f64,

    /// Also search reverse complement strand for primer matches
    #[arg(long, default_value_t = false)]
    pub search_rc: bool,

    /// Number of threads to use (0 = auto-detect)
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    // ==================== POOL OPTIONS ====================
    /// Treat primers as a pool of individual primers (all-vs-all matching).
    /// When set, -p format is "name:sequence" instead of "name:forward:reverse".
    #[arg(long, default_value_t = false)]
    pub pool: bool,

    /// In pool mode, allow the same primer to match as both forward and reverse
    #[arg(long, default_value_t = false)]
    pub pool_self_match: bool,

    // ==================== PRODUCT OPTIONS ====================
    /// Minimum PCR product length (including primers unless --trim-primers)
    #[arg(long, default_value_t = 50)]
    pub min_len: usize,

    /// Maximum PCR product length (including primers unless --trim-primers)
    #[arg(long, default_value_t = 5000)]
    pub max_len: usize,

    /// Remove primer sequences from output (output region between primers only)
    #[arg(long, default_value_t = false)]
    pub trim_primers: bool,

    /// Maximum fraction of N bases allowed in product sequence (0.0 - 1.0)
    #[arg(long, default_value_t = 0.1, value_parser = clap::value_parser!(f64))]
    pub max_n_fraction: f64,

    // ==================== OUTPUT OPTIONS ====================
    /// Output FASTA file (use .gz extension for gzip compression)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Optional TSV output file with detailed statistics
    #[arg(long)]
    pub tsv: Option<PathBuf>,

    /// Verbosity level (-v, -vv, -vvv)
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Suppress progress bar output
    #[arg(short, long, default_value_t = false)]
    pub quiet: bool,

    /// Remove duplicate product sequences per reference
    #[arg(long, default_value_t = false)]
    pub remove_duplicates: bool,
}

impl Args {
    /// Parse command line arguments
    #[must_use]
    pub fn parse_args() -> Self {
        Args::parse()
    }

    /// Get the effective number of threads
    #[must_use]
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            num_cpus::get()
        } else {
            self.threads
        }
    }

    /// Check if output should be gzip compressed
    #[must_use]
    pub fn output_is_gzipped(&self) -> bool {
        self.output
            .as_ref()
            .and_then(|p| p.extension())
            .is_some_and(|ext| ext == "gz")
    }

    /// Check if progress bar should be shown
    #[must_use]
    pub fn show_progress(&self) -> bool {
        !self.quiet && self.verbose == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_values() {
        let args = Args::parse_from(["amplirust", "-i", "test.fasta", "-p", "test:ACGT:TGCA"]);
        assert_eq!(args.max_errors, 2);
        assert!((args.min_identity - 0.0).abs() < 0.001);
        assert_eq!(args.min_len, 50);
        assert_eq!(args.max_len, 5000);
        assert!(!args.circular);
        assert!(!args.search_rc);
        assert!(!args.trim_primers);
    }

    #[test]
    fn test_multiple_inputs() {
        let args = Args::parse_from([
            "amplirust",
            "-i",
            "a.fasta,b.fasta,c.fasta",
            "-p",
            "test:ACGT:TGCA",
        ]);
        assert_eq!(args.input.len(), 3);
    }

    #[test]
    fn test_effective_threads_auto() {
        let args = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA"]);
        assert_eq!(args.threads, 0);
        assert!(args.effective_threads() > 0);
    }

    #[test]
    fn test_effective_threads_explicit() {
        let args = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA", "-t", "4"]);
        assert_eq!(args.effective_threads(), 4);
    }

    #[test]
    fn test_output_is_gzipped_true() {
        let args = Args::parse_from([
            "amplirust",
            "-i",
            "t.fa",
            "-p",
            "p:ACGT:TGCA",
            "-o",
            "out.fasta.gz",
        ]);
        assert!(args.output_is_gzipped());
    }

    #[test]
    fn test_output_is_gzipped_false() {
        let args = Args::parse_from([
            "amplirust",
            "-i",
            "t.fa",
            "-p",
            "p:ACGT:TGCA",
            "-o",
            "out.fasta",
        ]);
        assert!(!args.output_is_gzipped());
    }

    #[test]
    fn test_output_is_gzipped_no_output() {
        let args = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA"]);
        assert!(!args.output_is_gzipped());
    }

    #[test]
    fn test_show_progress() {
        let args = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA"]);
        assert!(args.show_progress());

        let quiet = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA", "-q"]);
        assert!(!quiet.show_progress());

        let verbose = Args::parse_from(["amplirust", "-i", "t.fa", "-p", "p:ACGT:TGCA", "-v"]);
        assert!(!verbose.show_progress());
    }
}
