//! Amplirust: In-silico PCR primer matching and product extraction tool
//!
//! This library provides functionality for:
//! - Parsing FASTA files (with gzip support)
//! - Approximate primer matching using IUPAC codes
//! - PCR product extraction with support for circular genomes
//! - Output in FASTA format with optional TSV statistics

pub mod cli;
pub mod errors;
pub mod genbank;
pub mod input;
pub mod matcher;
pub mod output;
pub mod pcr;
pub mod primer;
pub mod utils;

// Re-export main types for convenience
pub use cli::Args;
pub use input::SequenceRecord;
pub use matcher::{MatchConfig, PrimerMatch, PrimerMatcher};
pub use output::RunSummary;
pub use pcr::{PcrConfig, PcrProduct};
pub use primer::PrimerPair;
