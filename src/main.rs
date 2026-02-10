use anyhow::{Context, Result};
use log::LevelFilter;

use amplirust::cli::Args;
use amplirust::input::{expand_input_patterns, read_sequences_from_file, streaming_reader};
use amplirust::matcher::MatchConfig;
use amplirust::output::{FastaWriter, RunSummary, TsvWriter, write_fasta_record};
use amplirust::pcr::{PcrConfig, canonical_sequence, find_pcr_products};
use amplirust::primer::parse_primers;
use std::io::Write;
use std::sync::Arc;

fn main() -> Result<()> {
    // Parse command line arguments
    let args = Args::parse_args();

    // Initialize logging
    init_logging(args.verbose);

    // Run the main application
    run(args)
}

fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp(None)
        .format_target(false)
        .init();
}

fn run(args: Args) -> Result<()> {
    let show_progress = args.show_progress();

    log::info!("Amplirust v{}", env!("CARGO_PKG_VERSION"));

    if !(0.0..=1.0).contains(&args.max_n_fraction) {
        anyhow::bail!(
            "Invalid --max-n-fraction value {} (expected 0.0 - 1.0)",
            args.max_n_fraction
        );
    }

    if args.min_len > args.max_len {
        anyhow::bail!(
            "--min-len ({}) cannot exceed --max-len ({})",
            args.min_len,
            args.max_len
        );
    }

    if args.min_identity > 0.0 {
        if !(0.0..=1.0).contains(&args.min_identity) {
            anyhow::bail!(
                "Invalid --min-identity value {} (expected 0.0 - 1.0)",
                args.min_identity
            );
        }
        eprintln!(
            "Warning: --min-identity {:.2} is set; matches below {:.1}% identity will be filtered",
            args.min_identity,
            args.min_identity * 100.0
        );
    }

    // Configure thread pool
    let threads = args.effective_threads();
    log::info!("Using {threads} threads");

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Parse primers
    log::info!("Parsing primers...");
    let primers = parse_primers(&args.primers).context("Failed to parse primers")?;
    log::info!("Loaded {} primer pair(s)", primers.len());

    for primer in &primers {
        log::debug!(
            "Primer '{}': fwd={} ({}bp), rev={} ({}bp)",
            primer.name,
            String::from_utf8_lossy(&primer.forward),
            primer.forward_len(),
            String::from_utf8_lossy(&primer.reverse),
            primer.reverse_len()
        );
    }

    // Expand input file patterns
    log::info!("Finding input files...");
    let input_files =
        expand_input_patterns(&args.input).context("Failed to expand input patterns")?;

    // Configure PCR product detection
    let match_config = MatchConfig {
        max_errors: args.max_errors,
        min_identity: args.min_identity,
        search_rc: args.search_rc,
    };

    let pcr_config = PcrConfig {
        match_config,
        min_len: args.min_len,
        max_len: args.max_len,
        circular: args.circular,
        trim_primers: args.trim_primers,
        max_n_fraction: args.max_n_fraction,
    };

    log::info!("Searching for PCR products...");
    log::debug!("  Max errors: {}", args.max_errors);
    if args.min_identity > 0.0 {
        log::debug!("  Min identity: {:.1}%", args.min_identity * 100.0);
    } else {
        log::debug!("  Min identity: disabled");
    }
    log::debug!("  Product length: {}-{} bp", args.min_len, args.max_len);
    log::debug!("  Circular mode: {}", args.circular);
    log::debug!("  Search RC: {}", args.search_rc);
    log::debug!("  Trim primers: {}", args.trim_primers);
    log::debug!("  Max N fraction: {:.2}", args.max_n_fraction);

    // Process files in parallel and stream products to a single writer
    use indicatif::{ProgressBar, ProgressStyle};
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::mpsc::sync_channel;

    #[derive(Debug)]
    struct WorkerBatch {
        products: Vec<amplirust::pcr::PcrProduct>,
        sequences: usize,
    }

    #[derive(Debug)]
    struct RunOutcome {
        summary: RunSummary,
        wrote_fasta: bool,
        wrote_tsv: bool,
    }

    let pb = if show_progress && input_files.len() > 1 {
        let pb = ProgressBar::new(input_files.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files ({eta})")
                .unwrap()
                .progress_chars("#>-"),
        );
        pb.set_message("Reading files");
        Some(pb)
    } else {
        None
    };

    let completed = AtomicUsize::new(0);
    let (tx, rx) = sync_channel::<WorkerBatch>(threads * 2);

    let primers = Arc::new(primers);
    let pcr_config = Arc::new(pcr_config);
    let output_path = args.output.clone();
    let tsv_path = args.tsv.clone();
    let remove_duplicates = args.remove_duplicates;
    let total_primers = primers.len();
    let consumer_threads = threads;

    let consumer_handle = std::thread::spawn(move || -> Result<RunOutcome> {
        use std::collections::{HashMap, HashSet};
        let mut fasta_writer: Option<FastaWriter> = None;
        let mut tsv_writer: Option<TsvWriter> = None;
        let mut stdout_writer = if output_path.is_none() {
            Some(std::io::BufWriter::new(std::io::stdout().lock()))
        } else {
            None
        };

        let mut primer_counts: HashMap<String, usize> = HashMap::new();
        let mut ref_counts: HashMap<String, usize> = HashMap::new();
        let mut total_sequences = 0usize;
        let mut total_products = 0usize;
        let mut seen: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();

        let mut wrote_fasta = false;
        let mut wrote_tsv = false;

        while let Ok(batch) = rx.recv() {
            total_sequences += batch.sequences;
            for mut product in batch.products {
                if remove_duplicates {
                    let canonical = canonical_sequence(&product.sequence);
                    let entry = seen.entry(product.reference_header.clone()).or_default();
                    if entry.contains(&canonical) {
                        continue;
                    }
                    entry.insert(canonical);
                }

                let counter = ref_counts
                    .entry(product.reference_header.clone())
                    .or_insert(0);
                *counter += 1;
                product.case_number = *counter;

                total_products += 1;
                *primer_counts
                    .entry(product.primer_name.clone())
                    .or_insert(0) += 1;

                if let Some(ref path) = output_path {
                    if fasta_writer.is_none() {
                        fasta_writer = Some(FastaWriter::new(path, consumer_threads)?);
                    }
                    if let Some(ref mut writer) = fasta_writer {
                        writer.write_product(&product)?;
                        wrote_fasta = true;
                    }
                } else if let Some(ref mut writer) = stdout_writer {
                    write_fasta_record(writer, &product)?;
                }

                if let Some(ref path) = tsv_path {
                    if tsv_writer.is_none() {
                        tsv_writer = Some(TsvWriter::new(path)?);
                    }
                    if let Some(ref mut writer) = tsv_writer {
                        writer.write_product(&product)?;
                        wrote_tsv = true;
                    }
                }
            }
        }

        if let Some(writer) = fasta_writer {
            writer.finish()?;
        }
        if let Some(mut writer) = stdout_writer {
            writer.flush()?;
        }
        if let Some(writer) = tsv_writer {
            writer.finish()?;
        }

        let summary = RunSummary::from_counts(
            total_sequences,
            total_primers,
            total_products,
            primer_counts,
            ref_counts,
        );

        Ok(RunOutcome {
            summary,
            wrote_fasta,
            wrote_tsv,
        })
    });

    // Determine parallelism strategy
    let use_sequence_parallelism = input_files.len() == 1 && threads > 1;

    log::info!("Reading sequences from {} file(s)...", input_files.len());

    let results: Vec<Result<()>> = if use_sequence_parallelism {
        // Single file: use streaming sequence-level parallelism
        // Reader thread parses lazily, workers process in parallel via bounded channel
        log::info!("Using streaming sequence-level parallelism for single file");

        let file = &input_files[0];

        // Bounded channel from reader to workers (~200 sequences max in flight)
        // This provides backpressure when workers are slow
        const SEQUENCE_BUFFER_SIZE: usize = 200;
        let (seq_tx, seq_rx) =
            crossbeam_channel::bounded::<amplirust::input::SequenceRecord>(SEQUENCE_BUFFER_SIZE);

        // Progress bar (spinner mode - total unknown until file fully read)
        let processed_count = Arc::new(AtomicUsize::new(0));
        let seq_pb = if show_progress {
            let pb = ProgressBar::new_spinner();
            pb.set_style(
                ProgressStyle::default_spinner()
                    .template("{spinner:.green} [{elapsed_precise}] {pos} sequences processed")
                    .unwrap(),
            );
            Some(Arc::new(pb))
        } else {
            None
        };

        // Reader thread: parses sequences lazily, sends to workers
        let file_clone = file.clone();
        let reader_handle = std::thread::spawn(move || -> Result<()> {
            let reader = streaming_reader(&file_clone)?;
            for result in reader {
                let record = result?;
                if seq_tx.send(record).is_err() {
                    break; // Workers done or channel closed
                }
            }
            Ok(())
        });

        // Worker threads: pull from seq_rx, process, send to output channel
        // Drop the original tx before scope so consumer knows when workers are done
        let worker_tx = tx.clone();
        drop(tx);

        rayon::scope(|s| {
            for _ in 0..threads {
                let seq_rx = seq_rx.clone();
                let worker_tx = worker_tx.clone();
                let primers = Arc::clone(&primers);
                let pcr_config = Arc::clone(&pcr_config);
                let seq_pb = seq_pb.clone();
                let processed_count = Arc::clone(&processed_count);

                s.spawn(move |_| {
                    while let Ok(record) = seq_rx.recv() {
                        let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
                        if let Some(ref pb) = seq_pb {
                            pb.set_position(count as u64);
                        }

                        let mut products = Vec::new();
                        for primer in primers.iter() {
                            products.extend(find_pcr_products(
                                &record,
                                primer,
                                pcr_config.as_ref(),
                            ));
                        }

                        // Send batch of 1 sequence's products (maintains streaming flow)
                        let _ = worker_tx.send(WorkerBatch {
                            products,
                            sequences: 1,
                        });
                    }
                });
            }
        });

        // Drop worker_tx to signal consumer that workers are done
        drop(worker_tx);

        // Wait for reader thread
        reader_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Reader thread panicked"))??;

        let total_sequences = processed_count.load(Ordering::Relaxed);
        if let Some(ref pb) = seq_pb {
            pb.finish_with_message(format!("{total_sequences} sequences processed"));
        }

        // tx already dropped, return empty results (no per-file errors in streaming mode)
        return match consumer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Writer thread panicked"))?
        {
            Ok(outcome) => {
                // Print summary and output info
                outcome.summary.print_stderr();
                outcome.summary.log_summary();

                if outcome.summary.total_products == 0 {
                    log::warn!("No PCR products found matching the criteria");
                    return Ok(());
                }

                if show_progress {
                    if let Some(ref output_path) = args.output
                        && outcome.wrote_fasta
                    {
                        eprintln!("Output written to: {}", output_path.display());
                    }
                    if let Some(ref tsv_path) = args.tsv
                        && outcome.wrote_tsv
                    {
                        eprintln!("TSV written to: {}", tsv_path.display());
                    }
                }

                log::info!("Done!");
                Ok(())
            }
            Err(e) => Err(e),
        };
    } else {
        // Multiple files: use file-level parallelism
        input_files
            .par_iter()
            .map(|file| {
                let tx = tx.clone();
                log::debug!("Reading sequences from: {}", file.display());
                let records = read_sequences_from_file(file)
                    .with_context(|| format!("Failed to read input file: {}", file.display()))?;

                let mut products = Vec::new();
                for record in &records {
                    for primer in primers.iter() {
                        products.extend(find_pcr_products(record, primer, pcr_config.as_ref()));
                    }
                }

                tx.send(WorkerBatch {
                    products,
                    sequences: records.len(),
                })
                .context("Failed to send products to writer")?;

                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                if let Some(ref pb) = pb {
                    pb.set_position(done as u64);
                }

                Ok(())
            })
            .collect()
    };

    drop(tx);

    if let Some(pb) = pb {
        pb.finish_with_message("Files processed");
    }

    let mut worker_error = None;
    for result in results {
        if let Err(err) = result {
            worker_error = Some(err);
            break;
        }
    }

    let outcome = consumer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    if let Some(err) = worker_error {
        return Err(err);
    }

    // Always print summary to stderr (unless quiet mode and no products)
    outcome.summary.print_stderr();
    outcome.summary.log_summary();

    if outcome.summary.total_products == 0 {
        log::warn!("No PCR products found matching the criteria");
        return Ok(());
    }

    if show_progress {
        if let Some(ref output_path) = args.output
            && outcome.wrote_fasta
        {
            eprintln!("Output written to: {}", output_path.display());
        }
        if let Some(ref tsv_path) = args.tsv
            && outcome.wrote_tsv
        {
            eprintln!("TSV written to: {}", tsv_path.display());
        }
    }

    log::info!("Done!");
    Ok(())
}
