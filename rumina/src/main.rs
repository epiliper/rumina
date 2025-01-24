#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, GroupingMethod};
use crate::group_report::GroupReport;
use crate::process::{gather_files, process};
use indexmap::IndexMap;
use std::fs::create_dir;
use std::path::Path;

use log::LevelFilter;
use rayon::ThreadPoolBuilder;

mod args;
mod bottomhash;
mod deduplicator;
mod group_report;
mod grouper;
mod main_dedup;
mod merge;
mod merge_report;
mod pair_merger;
mod process;
mod progbars;
mod read_picker;
mod readkey;
mod realign;
mod utils;
mod window_processor;

fn main() {
    let args = parse_args();
    let input_file = &args.input;

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("ERROR: Invalid number of threads specified!");

    simple_logging::log_to_file("rumina_group.log", LevelFilter::Info)
        .expect("Failed to create log for grouping command!");

    if args.merge_pairs.is_some() {
        simple_logging::log_to_file("rumina_merge.log", LevelFilter::Info)
            .expect("Failed to create log for merge command!");
    }

    if !Path::exists(Path::new(&args.outdir)) {
        create_dir(&args.outdir).unwrap();
    }

    let infiles = gather_files(input_file);

    for file in infiles {
        process(file, &args);
    }
}
