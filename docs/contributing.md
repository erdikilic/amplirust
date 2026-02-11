# Contributing

Contributions to Amplirust are welcome!

## Development Setup

```bash
# Clone the repository
git clone https://github.com/erdikilic/amplirust.git
cd amplirust

# Build
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run tests
cargo test
```

## Code Quality

CI enforces a quality gate on every push and pull request. Before submitting changes, ensure your code passes all checks locally:

```bash
# Format code
cargo fmt --all

# Run linter (pedantic lints are configured in Cargo.toml)
cargo clippy --all-targets -- -D warnings

# Audit dependencies for security vulnerabilities and license issues
cargo deny check

# Run tests with both parser backends
cargo test
cargo test --no-default-features --features parser_needletail
```

For performance-sensitive changes, run the Criterion benchmarks to check for regressions:

```bash
cargo bench
```

## Feature Flags

Amplirust uses feature flags for FASTA parser selection:

- `parser_seqio` (default) -- uses the [seq_io](https://github.com/markschl/seq_io) parser
- `parser_needletail` -- uses the [needletail](https://github.com/onecodex/needletail) parser

Test with both parsers when making changes to input handling:

```bash
cargo test
cargo test --no-default-features --features parser_needletail
```

## Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-change`)
3. Make your changes and add tests
4. Ensure all CI checks pass: `cargo fmt --all`, `cargo clippy --all-targets -- -D warnings`, `cargo deny check`, and `cargo test` (with both parser features)
5. Open a pull request against `main`

## License

Amplirust is licensed under the [MIT License](https://github.com/erdikilic/amplirust/blob/main/LICENSE).
