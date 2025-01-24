#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, GroupingMethod};
use crate::group_report::GroupReport;
use crate::process::process;
use indexmap::IndexMap;
use std::fs::{create_dir, read_dir};
use std::path::Path;
use utils::get_file_ext;

use log::{error, LevelFilter};
use rayon::ThreadPoolBuilder;

use std::collections::HashMap;

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

    let mut infiles: HashMap<String, String> = HashMap::new();
    let inpath = Path::new(&input_file);
    if inpath.is_dir() {
        for file in read_dir(inpath)
            .unwrap()
            .filter_map(|f| match f {
                Err(e) => {
                    error!("Unable to read file: {e}. Won't be used in processing.");
                    None
                }
                Ok(f) => Some(f),
            })
            .filter(|f| !f.path().is_dir() && get_file_ext(f) == "bam")
        {
            infiles.insert(
                file.path().to_str().unwrap().to_string(),
                file.file_name().to_str().unwrap().to_string(),
            );
        }
    } else {
        infiles.insert(
            inpath.to_str().unwrap().to_string(),
            inpath.file_name().unwrap().to_str().unwrap().to_string(),
        );
    }

    if !Path::exists(Path::new(&args.outdir)) {
        create_dir(&args.outdir).unwrap();
    }

    for file in infiles {
        process(file, &args);
    }
}
