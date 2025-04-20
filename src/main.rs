#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, GroupingMethod};
use crate::cli::*;
use crate::group_report::GroupReport;
use crate::io::{gather_files, process_all};
use indexmap::IndexMap;
use std::fs::create_dir;
use std::path::Path;

use anyhow::{Context, Error};
use log::LevelFilter;
use rayon::ThreadPoolBuilder;

mod args;
mod cli;
mod deduplicator;
mod group;
mod group_report;
mod grouper;
mod io;
mod merge;
mod merge_report;
mod ngram;
mod pair_merger;
mod process;
mod processor;
mod progbars;
mod read_picker;
mod read_store;
mod readkey;
mod realign;
mod record;
mod utils;

fn main() -> Result<(), Error> {
    let args = parse_args();
    let input_file = &args.input;

    print_logo();
    print_init(&args);

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .with_context(|| "Thread pool building failed")?;

    simple_logging::log_to_file("rumina_group.log", LevelFilter::Info)?;

    if args.merge_pairs.is_some() {
        simple_logging::log_to_file("rumina_merge.log", LevelFilter::Info)?;
    }

    if !Path::exists(Path::new(&args.outdir)) {
        create_dir(&args.outdir)
            .with_context(|| format!("Unable to create output directory {}", &args.outdir))?;
    }

    let infiles = gather_files(input_file)?;

    process_all(&args, infiles)
        .into_iter()
        .map(|r| r.err())
        .flatten()
        .for_each(|e| eprintln! {"{:?}\n--", e});

    Ok(())
}
