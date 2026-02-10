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
        /// File path that failed the writability check.
        path: PathBuf,
        /// Underlying I/O error from the write probe.
        #[source]
        source: std::io::Error,
    },

    /// Parent directory of the output path does not exist.
    #[error("output directory '{}' does not exist", path.display())]
    OutputDirMissing {
        /// Output file path whose parent directory is missing.
        path: PathBuf,
    },

    /// Decompressed input exceeded the configured size limit.
    #[error(
        "decompressed size of '{}' exceeds limit of {limit} bytes",
        path.display()
    )]
    DecompressionLimitExceeded {
        /// Path to the compressed input file.
        path: PathBuf,
        /// Maximum allowed decompressed size in bytes.
        limit: u64,
    },

    /// CSV/TSV primer file has a formatting issue.
    #[error("invalid CSV format in '{}': {detail}", path.display())]
    CsvFormat {
        /// Path to the malformed CSV/TSV file.
        path: PathBuf,
        /// Human-readable description of the formatting problem.
        detail: String,
    },

    /// A single line in an input file exceeds the safety limit.
    #[error(
        "line length {len} in '{}' exceeds limit of {limit} bytes",
        path.display()
    )]
    LineTooLong {
        /// Path to the file containing the oversized line.
        path: PathBuf,
        /// Actual length of the offending line in bytes.
        len: usize,
        /// Configured maximum line length in bytes.
        limit: usize,
    },
}
