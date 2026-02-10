use std::hint::black_box;
use std::path::PathBuf;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use amplirust::input::SequenceRecord;
use amplirust::matcher::MatchConfig;
use amplirust::pcr::{PcrConfig, find_pcr_products};
use amplirust::primer::PrimerPair;
use amplirust::utils::reverse_complement;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// Generate a deterministic random DNA sequence of `len` bytes using the given seed.
fn generate_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len).map(|_| bases[rng.random_range(0..4usize)]).collect()
}

/// Generate a sequence with forward and reverse primers planted at known positions.
///
/// Plants the forward primer at position `len/4` (or 1000 if `len >= 4000`)
/// and the reverse primer RC at position `3*len/4` (or `len-1000` if `len >= 4000`).
/// This ensures at least one PCR product will be found.
fn generate_sequence_with_primers(len: usize, fwd: &[u8], rev_rc: &[u8], seed: u64) -> Vec<u8> {
    let mut seq = generate_dna(len, seed);

    let fwd_pos = if len >= 4000 { 1000 } else { len / 4 };
    let rev_pos = if len >= 4000 { len - 1000 } else { 3 * len / 4 };

    // Plant forward primer
    let fwd_end = (fwd_pos + fwd.len()).min(len);
    seq[fwd_pos..fwd_end].copy_from_slice(&fwd[..fwd_end - fwd_pos]);

    // Plant reverse primer RC
    let rev_end = (rev_pos + rev_rc.len()).min(len);
    seq[rev_pos..rev_end].copy_from_slice(&rev_rc[..rev_end - rev_pos]);

    seq
}

fn bench_full_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_pipeline");

    // 16S universal primers
    let fwd: &[u8] = b"AGAGTTTGATCMTGGCTCAG";
    let rev: &[u8] = b"TACGGYTACCTTGTTACGACTT";
    let rev_rc = reverse_complement(rev);

    let primer = PrimerPair::new("16S", fwd, rev).unwrap();

    let config = PcrConfig {
        match_config: MatchConfig {
            max_errors: 2,
            min_identity: 0.8,
            search_rc: true,
        },
        min_len: 50,
        max_len: 5000,
        circular: false,
        trim_primers: false,
        max_n_fraction: 1.0,
    };

    let sizes: &[usize] = &[50_000, 200_000, 1_000_000];

    for &size in sizes {
        let seq = generate_sequence_with_primers(size, fwd, &rev_rc, 42);
        let record = SequenceRecord {
            header: format!("bench_seq_{size}"),
            sequence: seq,
            source_file: PathBuf::from("bench.fasta"),
        };

        group.throughput(Throughput::Bytes(size as u64));

        group.bench_with_input(BenchmarkId::from_parameter(size), &record, |b, record| {
            b.iter(|| find_pcr_products(black_box(record), black_box(&primer), black_box(&config)));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_full_pipeline);
criterion_main!(benches);
