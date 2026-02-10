#![no_main]

use libfuzzer_sys::fuzz_target;
use seq_io::fasta::{Reader, Record};

fuzz_target!(|data: &[u8]| {
    let mut reader = Reader::new(data);
    while let Some(result) = reader.next() {
        match result {
            Ok(record) => {
                let _ = record.id();
                let _ = record.full_seq();
            }
            Err(_) => break,
        }
    }
});
