#![no_main]

use libfuzzer_sys::fuzz_target;

use amplirust::genbank::parse_genbank_slice;

fuzz_target!(|data: &[u8]| {
    let records = parse_genbank_slice(data);
    for rec in &records {
        let _ = &rec.name;
        let _ = &rec.accession;
        let _ = &rec.definition;
        let _ = rec.is_circular;
        let _ = rec.seq.len();
    }
});
