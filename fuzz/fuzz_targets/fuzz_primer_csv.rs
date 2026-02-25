#![no_main]

use libfuzzer_sys::fuzz_target;
use std::io::Write;
use tempfile::NamedTempFile;

fuzz_target!(|data: &[u8]| {
    let Ok(mut temp) = NamedTempFile::with_suffix(".csv") else {
        return;
    };
    if temp.write_all(data).is_err() {
        return;
    }
    if temp.flush().is_err() {
        return;
    }
    let _ = amplirust::primer::parse_primers(&temp.path().to_string_lossy());
});
