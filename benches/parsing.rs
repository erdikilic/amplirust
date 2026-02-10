// Placeholder: FASTA/GenBank parsing benchmarks (Plan 06-02)
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_parsing(_c: &mut Criterion) {
    // TODO: Implement in Plan 06-02
}

criterion_group!(benches, bench_parsing);
criterion_main!(benches);
