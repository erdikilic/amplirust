#![warn(missing_docs)]
//! Amplirust: In-silico PCR primer matching and product extraction tool
//!
//! This library provides functionality for:
//! - Parsing FASTA files (with gzip support)
//! - Approximate primer matching using IUPAC codes
//! - PCR product extraction with support for circular genomes
//! - Output in FASTA format with optional TSV statistics

/// Command-line argument definitions and parsing.
pub mod cli;
/// Structured error types for input validation and runtime failures.
pub mod errors;
/// GenBank flat-file format parser with streaming and slice-based modes.
pub mod genbank;
/// Sequence input handling: FASTA/GenBank detection, gzip decompression, file expansion.
pub mod input;
/// Approximate primer matching using IUPAC-aware edit distance.
pub mod matcher;
/// Output formatting for FASTA sequences and TSV statistics.
pub mod output;
/// In-silico PCR product detection with forward, reverse-complement, and circular support.
pub mod pcr;
/// Multi-threaded pipeline orchestrating file parsing, matching, and output.
pub mod pipeline;
/// Primer pair parsing from CSV files and inline strings.
pub mod primer;
/// DNA sequence utilities: complement, reverse complement, circular extension.
pub mod utils;

// Re-export main types for convenience
pub use cli::Args;
pub use input::SequenceRecord;
pub use matcher::{MatchConfig, PrimerMatch, PrimerMatcher};
pub use output::RunSummary;
pub use pcr::{PcrConfig, PcrProduct, find_pool_products};
pub use primer::{Primer, PrimerPair};
