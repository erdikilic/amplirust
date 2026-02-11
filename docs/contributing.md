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

## Testing

### Unit & integration tests

```bash
cargo test
cargo test --no-default-features --features parser_needletail
```

### Fuzz testing

Requires a nightly toolchain and [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz):

```bash
cargo install cargo-fuzz

# Run a target (runs until stopped with Ctrl-C)
cargo +nightly fuzz run fuzz_fasta

# Available targets: fuzz_fasta, fuzz_genbank_streaming, fuzz_genbank_slice,
#                    fuzz_primer_csv, fuzz_primer_matcher
cargo +nightly fuzz list
```

### Mutation testing

Requires [cargo-mutants](https://github.com/sourcefrog/cargo-mutants):

```bash
cargo install cargo-mutants

cargo mutants
```

Configuration is in `.cargo/mutants.toml`.

### Benchmarks

Requires [Criterion](https://bheisler.github.io/criterion.rs/book/) (already a dev-dependency):

```bash
# Run all benchmarks
cargo bench

# Run a specific suite
cargo bench --bench matching
cargo bench --bench parsing
cargo bench --bench pipeline
cargo bench --bench circular
```

### Validation against UCSC isPCR

The script `scripts/validate_ispcr.sh` compares amplirust output against the gold-standard UCSC isPCR tool. It requires `isPcr` on your PATH and downloads the E. coli K-12 genome:

```bash
bash scripts/validate_ispcr.sh
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
