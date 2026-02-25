use anyhow::Result;
use log::LevelFilter;

use amplirust::cli::Args;

fn main() -> Result<()> {
    // Parse command line arguments
    let args = Args::parse_args();

    // Initialize logging
    init_logging(args.verbose);

    // Run the main application
    amplirust::pipeline::run(&args)
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
