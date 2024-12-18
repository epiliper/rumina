#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::main_dedup::{init_processor, process_chunks};
use crate::report::{BarcodeTracker, GroupReport};
use clap::ValueEnum;
use colored::Colorize;
use indexmap::IndexMap;
use std::hash::{DefaultHasher, Hash, Hasher};

use clap::Parser;
use rayon::ThreadPoolBuilder;

use parking_lot::Mutex;
use std::sync::Arc;

mod bottomhash;
mod deduplicator;
mod grouper;
mod main_dedup;
mod progbars;
mod read_picker;
mod readkey;
mod report;
mod utils;
mod window_processor;

#[derive(ValueEnum, Debug, Clone)]
enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}

#[derive(Parser, Debug)]
#[command(term_width = 0)]
struct Args {
    #[arg(long = "in", index = 1)]
    input: String,

    #[arg(long = "out", index = 2)]
    output: String,

    #[arg(long = "sep", index = 3)]
    separator: String,

    #[arg(long = "group", index = 4)]
    grouping_method: GroupingMethod,

    #[arg(long = "threads", index = 5)]
    threads: usize,

    #[arg(long = "split_window", index = 6)]
    split_window: Option<i64>,

    #[arg(long = "length")]
    length: bool,

    #[arg(long = "only-group")]
    only_group: bool,

    #[arg(long = "singletons")]
    singletons: bool,

    #[arg(long = "track-umis")]
    track_barcodes: bool,

    #[arg(long = "r1-only")]
    r1_only: bool,
}

fn main() {
    let args = Args::parse();
    let input_file = args.input;
    let output_file = &args.output;
    let separator = args.separator;
    let grouping_method = args.grouping_method;
    let split_window = args.split_window;
    let length = args.length;
    let only_group = args.only_group;
    let mut track_barcodes = None;
    let singletons = args.singletons;

    if args.track_barcodes {
        track_barcodes = Some(output_file)
    }

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        // .build()
        .build_global()
        .expect("ERROR: Invalid number of threads specified!");

    // create tag seed based on input file name
    let mut hasher = DefaultHasher::new();
    input_file.hash(&mut hasher);
    let seed = hasher.finish();

    // create deduplication report
    let min_maxes: Arc<Mutex<GroupReport>> = Arc::new(Mutex::new(GroupReport::new()));
    let barcode_tracker: Arc<Mutex<BarcodeTracker>> =
        Arc::new(Mutex::new(BarcodeTracker::new(&output_file)));

    let (bam_reader, bam_writer, mut read_handler) = init_processor(
        input_file,
        output_file.to_string(),
        grouping_method,
        args.threads,
        split_window,
        length,
        only_group,
        singletons,
        track_barcodes,
        args.r1_only,
        min_maxes.clone(),
        barcode_tracker,
        seed,
    );

    // holds filtered reads awaiting writing to output bam file
    // do grouping and processing
    process_chunks(&mut read_handler, bam_reader, &separator, bam_writer);
    let num_reads_in = read_handler.read_counter;

    drop(read_handler);

    // do final report
    let mut group_report = Arc::try_unwrap(min_maxes).unwrap().into_inner();
    group_report.num_reads_input_file = num_reads_in;

    // report on min and max number of reads per group
    // this creates minmax.txt
    if group_report.min_reads != i64::MAX {
        println!("{}", "DONE".green());

        group_report.write_to_report_file(&output_file);
        println!("{}", group_report);
    }
}
