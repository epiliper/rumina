#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, GroupingMethod};
use crate::cli::*;
use crate::group_report::GroupReport;
use crate::process::{gather_files, FileProcess};
use indexmap::IndexMap;
use std::fs::create_dir;
use std::path::Path;

use log::LevelFilter;
use rayon::ThreadPoolBuilder;

mod args;
mod bam_io;
mod cli;
mod deduplicator;
mod group;
mod group_report;
mod grouper;
mod main_dedup;
mod merge;
mod merge_report;
mod ngram;
mod pair_merger;
mod process;
mod progbars;
mod read_picker;
mod read_store;
mod readkey;
mod realign;
mod utils;
mod window_processor;

fn main() {
    let args = parse_args();
    let input_file = &args.input;

    print_logo();
    print_init(&args);

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

    let num_files = infiles.len();
    for (i, (path, name)) in infiles.into_iter().enumerate() {
        print_file_info(&name, i + 1, num_files);
        process::BamFileProcess::init_from_args(&args, &path, &name).process();
    }
}
