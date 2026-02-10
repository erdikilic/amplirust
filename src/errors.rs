//! Structured error types for input validation and parse warnings.

use std::path::PathBuf;

use thiserror::Error;

/// Validation errors for input files and configuration.
///
/// Each variant carries enough context (file path, counts, limits) to produce
/// an actionable error message without requiring the caller to add extra context.
#[derive(Debug, Error)]
pub enum ValidationError {
    /// Output path exists but is not writable, or creation failed.
    #[error("output file '{}' is not writable: {source}", path.display())]
    OutputNotWritable {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// Parent directory of the output path does not exist.
    #[error("output directory '{}' does not exist", path.display())]
    OutputDirMissing { path: PathBuf },

    /// Decompressed input exceeded the configured size limit.
    #[error(
        "decompressed size of '{}' exceeds limit of {limit} bytes",
        path.display()
    )]
    DecompressionLimitExceeded { path: PathBuf, limit: u64 },

    /// CSV/TSV primer file has a formatting issue.
    #[error("invalid CSV format in '{}': {detail}", path.display())]
    CsvFormat { path: PathBuf, detail: String },

    /// A single line in an input file exceeds the safety limit.
    #[error(
        "line length {len} in '{}' exceeds limit of {limit} bytes",
        path.display()
    )]
    LineTooLong {
        path: PathBuf,
        len: usize,
        limit: usize,
    },
}
