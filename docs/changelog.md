# Changelog

## Unreleased

### Added

- `--max-decompression-size` CLI option to set a safety limit against decompression bombs (default 4 GiB, 0 = unlimited)
- Structured error handling with actionable messages for common issues (unwritable output, missing directories, decompression limits, CSV format errors)
- 260 tests: unit, integration, property-based (proptest), snapshot (insta), fuzz regression, and validation
- Criterion benchmarks for matching, parsing, pipeline, and circular genome operations
- Complete public API documentation enforced by `missing_docs` lint
- CI quality gate on every push and pull request (rustfmt, clippy pedantic, cargo-deny, coverage with Codecov)

### Changed

- Internal: extracted pipeline logic to `src/pipeline.rs`, unified GenBank parser, adaptive buffer sizing

## v0.2.0

Current release.

- SIMD-accelerated approximate primer matching via [sassy](https://github.com/RagnarGrootKoerkamp/sassy)
- Full IUPAC ambiguity code support in primers
- FASTA and GenBank input with automatic format detection
- Circular genome support with wrap-around product detection
- Reverse complement strand search
- Multi-threaded file reading, searching, and compression
- Gzip/BGZF transparent compression for input and output
- Flexible primer input via command line or CSV file
- TSV statistics output with alignment details (CIGAR strings, identity scores)
- Configurable product length filtering and N-base fraction limits
- Duplicate product removal (canonicalized by reverse complement)
- Progress bar with run summary
- Two FASTA parser backends: seq_io (default) and needletail
- Available via Bioconda for Linux and macOS (x86_64, aarch64)
