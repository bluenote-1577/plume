use crate::cmdline::*;
use log::*;
use rayon::prelude::*;

pub fn quant(args: QuantArgs){

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else {
        level = log::LevelFilter::Info;
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    if args.files.is_empty()
    {
        error!("No input sequences found. Exiting.");
        std::process::exit(1);
    }
}

