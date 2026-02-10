#![no_main]

use libfuzzer_sys::fuzz_target;
use std::io::Cursor;

use amplirust::genbank::GenbankReader;

fuzz_target!(|data: &[u8]| {
    let reader = GenbankReader::new(Cursor::new(data));
    for result in reader {
        let _ = result;
    }
});
