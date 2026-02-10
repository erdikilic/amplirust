use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use amplirust::matcher::{MatchConfig, PrimerMatcher};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// Generate a deterministic random DNA sequence of `len` bytes using the given seed.
fn generate_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len).map(|_| bases[rng.random_range(0..4usize)]).collect()
}

fn bench_primer_matching(c: &mut Criterion) {
    let mut group = c.benchmark_group("primer_matching");

    // 16S universal forward primer (27F): 20 bp with IUPAC ambiguity (M = A or C)
    let primer: &[u8] = b"AGAGTTTGATCMTGGCTCAG";

    let sizes: &[usize] = &[10_000, 100_000, 1_000_000];

    for &size in sizes {
        let seq = generate_dna(size, 42);

        group.throughput(Throughput::Bytes(size as u64));

        // Exact matching (no errors allowed)
        let exact_config = MatchConfig {
            max_errors: 0,
            min_identity: 1.0,
            search_rc: false,
        };

        group.bench_with_input(
            BenchmarkId::new("exact", size),
            &size,
            |b, _| {
                let mut matcher = PrimerMatcher::new(exact_config.clone());
                b.iter(|| matcher.find_matches(black_box(primer), black_box(&seq)));
            },
        );

        // Approximate matching (up to 2 errors, 80% identity)
        let approx_config = MatchConfig {
            max_errors: 2,
            min_identity: 0.8,
            search_rc: false,
        };

        group.bench_with_input(
            BenchmarkId::new("approx_2err", size),
            &size,
            |b, _| {
                let mut matcher = PrimerMatcher::new(approx_config.clone());
                b.iter(|| matcher.find_matches(black_box(primer), black_box(&seq)));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_primer_matching);
criterion_main!(benches);
