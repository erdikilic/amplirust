use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use amplirust::genbank::parse_genbank_slice;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// Generate a deterministic random DNA sequence of `len` bytes using the given seed.
fn generate_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len).map(|_| bases[rng.random_range(0..4usize)]).collect()
}

/// Generate synthetic `GenBank` data with `num_records` records, each with `seq_len` bp.
fn generate_genbank_data(num_records: usize, seq_len: usize, seed: u64) -> Vec<u8> {
    let mut data = Vec::new();

    for i in 0..num_records {
        // LOCUS line
        data.extend_from_slice(
            format!("LOCUS       SEQ{i:<13} {seq_len} bp    DNA     linear   UNK\n").as_bytes(),
        );
        // DEFINITION line
        data.extend_from_slice(format!("DEFINITION  Benchmark sequence {i}\n").as_bytes());
        // ACCESSION line
        data.extend_from_slice(format!("ACCESSION   SEQ{i}\n").as_bytes());
        // ORIGIN header
        data.extend_from_slice(b"ORIGIN\n");

        // Generate the sequence
        let seq = generate_dna(seq_len, seed + i as u64);

        // Write in GenBank ORIGIN format: line number (9-char right-aligned),
        // then up to 6 groups of 10 bases separated by spaces, 60 bases per line
        let mut pos = 0;
        while pos < seq.len() {
            // Line number (1-based, 9 chars right-aligned)
            let line_num = pos + 1;
            data.extend_from_slice(format!("{line_num:>9}").as_bytes());

            // Up to 60 bases per line, in groups of 10
            let line_end = (pos + 60).min(seq.len());
            let mut col = pos;
            while col < line_end {
                data.push(b' ');
                let group_end = (col + 10).min(line_end);
                data.extend_from_slice(&seq[col..group_end]);
                col = group_end;
            }
            data.push(b'\n');
            pos = line_end;
        }

        // Record terminator
        data.extend_from_slice(b"//\n");
    }

    data
}

fn bench_genbank_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("genbank_parsing");

    let cases: &[(&str, usize, usize)] = &[
        ("1x100kb", 1, 100_000),
        ("100x1kb", 100, 1_000),
        ("10x10kb", 10, 10_000),
    ];

    for &(name, num_records, seq_len) in cases {
        let data = generate_genbank_data(num_records, seq_len, 42);

        group.throughput(Throughput::Bytes(data.len() as u64));

        group.bench_with_input(BenchmarkId::from_parameter(name), &data, |b, data| {
            b.iter(|| parse_genbank_slice(black_box(data)));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_genbank_parsing);
criterion_main!(benches);
