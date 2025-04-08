#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::cli::args::{parse_args, GroupingMethod};
use crate::cluster::group_report::GroupReport;
use crate::process::{gather_files, process};
use indexmap::IndexMap;
use std::fs::create_dir;
use std::path::Path;

use log::LevelFilter;
use rayon::ThreadPoolBuilder;

mod bam_io;
mod bottomhash;
mod cli;
mod cluster;
mod group_handler;
mod main_dedup;
mod pair_merge;
mod process;
mod read_picker;
mod readkey;
mod utils;
mod window_processor;

fn main() {
    let args = parse_args();
    let input_file = &args.input;

    cli::info::print_logo();
    cli::info::print_init(&args);

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
    for (i, file) in infiles.into_iter().enumerate() {
        cli::info::print_file_info(&file.0, i + 1, num_files);
        process(file, &args);
    }
}
