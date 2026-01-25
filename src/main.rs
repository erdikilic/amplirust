use anyhow::{Context, Result};
use log::LevelFilter;

use amplirust::cli::Args;
use amplirust::input::{expand_input_patterns, read_all_sequences};
use amplirust::matcher::MatchConfig;
use amplirust::output::{write_fasta, write_fasta_stdout, write_tsv, RunSummary};
use amplirust::pcr::{find_all_products, remove_duplicate_products_by_reference, PcrConfig};
use amplirust::primer::parse_primers;

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

    // Configure thread pool
    let threads = args.effective_threads();
    log::info!("Using {} threads", threads);
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Parse primers
    log::info!("Parsing primers...");
    let primers = parse_primers(&args.primers)
        .context("Failed to parse primers")?;
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
    let input_files = expand_input_patterns(&args.input)
        .context("Failed to expand input patterns")?;

    // Read all sequences (with progress bar)
    log::info!("Reading sequences from {} file(s)...", input_files.len());
    let sequences = read_all_sequences(&input_files, show_progress)
        .context("Failed to read input sequences")?;
    
    if sequences.is_empty() {
        log::warn!("No sequences found in input files");
        eprintln!("Warning: No sequences found in input files");
        return Ok(());
    }
    log::info!("Loaded {} sequence(s)", sequences.len());

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
    log::debug!("  Min identity: {:.1}%", args.min_identity * 100.0);
    log::debug!("  Product length: {}-{} bp", args.min_len, args.max_len);
    log::debug!("  Circular mode: {}", args.circular);
    log::debug!("  Search RC: {}", args.search_rc);
    log::debug!("  Trim primers: {}", args.trim_primers);
    log::debug!("  Max N fraction: {:.2}", args.max_n_fraction);

    // Find all PCR products (with progress bar)
    let mut products = find_all_products(&sequences, &primers, &pcr_config, show_progress);
    if args.remove_duplicates {
        let before = products.len();
        products = remove_duplicate_products_by_reference(products);
        let removed = before.saturating_sub(products.len());
        if removed > 0 {
            log::info!("Removed {} duplicate product(s) per reference", removed);
        }
    }

    // Generate summary
    let summary = RunSummary::from_products(&products, sequences.len(), primers.len());
    
    // Always print summary to stderr (unless quiet mode and no products)
    summary.print_stderr();
    
    // Also log for verbose mode
    summary.log_summary();

    if products.is_empty() {
        log::warn!("No PCR products found matching the criteria");
        return Ok(());
    }

    // Write output
    if let Some(ref output_path) = args.output {
        write_fasta(&products, output_path, threads)
            .with_context(|| format!("Failed to write output to {}", output_path.display()))?;
        if show_progress {
            eprintln!("Output written to: {}", output_path.display());
        }
    } else {
        // Write to stdout if no output file specified
        write_fasta_stdout(&products)
            .context("Failed to write output to stdout")?;
    }

    // Write TSV if requested
    if let Some(ref tsv_path) = args.tsv {
        write_tsv(&products, tsv_path)
            .with_context(|| format!("Failed to write TSV to {}", tsv_path.display()))?;
        if show_progress {
            eprintln!("TSV written to: {}", tsv_path.display());
        }
    }

    log::info!("Done!");
    Ok(())
}
