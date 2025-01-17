#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::main_dedup::{init_processor, process_chunks};
use crate::pair_merger::PairMerger;
use crate::report::GroupReport;
use colored::Colorize;
use indexmap::IndexMap;
use std::hash::{DefaultHasher, Hash, Hasher};

use clap::{Parser, Subcommand, ValueEnum};
use rayon::ThreadPoolBuilder;

use parking_lot::Mutex;
use std::sync::Arc;

mod bottomhash;
mod deduplicator;
mod grouper;
mod main_dedup;
mod merge;
mod merge_report;
mod pair_merger;
mod progbars;
mod read_picker;
mod readkey;
mod realign;
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
    #[command(subcommand)]
    command: Command,

    #[arg(long = "in", index = 1)]
    input: String,

    #[arg(long = "out", index = 2)]
    output: String,

    #[arg(long = "threads", index = 3)]
    threads: usize,

    #[arg(long = "split_window")]
    split_window: Option<i64>,
}

#[derive(Subcommand, Debug)]
enum Command {
    Merge {
        #[arg(long = "merge_pairs")]
        ref_fasta: String,

        #[arg(long = "min_overlap_bp")]
        min_overlap_bp: usize,
    },

    Group {
        #[arg(short = 's')]
        separator: String,

        #[arg(short = 'g')]
        grouping_method: GroupingMethod,

        #[arg(long = "length")]
        length: bool,

        #[arg(long = "only-group")]
        only_group: bool,

        #[arg(long = "singletons")]
        singletons: bool,

        #[arg(long = "r1-only")]
        r1_only: bool,
    },
}

fn main() {
    let args = Args::parse();
    let input_file = args.input;
    let output_file = &args.output;

    match args.command {
        Command::Group {
            separator,
            grouping_method,
            length,
            only_group,
            singletons,
            r1_only,
        } => {
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

            let (bam_reader, bam_writer, mut read_handler) = init_processor(
                input_file,
                output_file.to_string(),
                grouping_method,
                args.threads,
                args.split_window,
                length,
                only_group,
                singletons,
                r1_only,
                min_maxes.clone(),
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

        Command::Merge {
            ref_fasta,
            min_overlap_bp,
        } => {
            let mut merger = PairMerger {
                ref_fasta,
                min_overlap_bp,
                threads: args.threads,
                infile: input_file.to_string(),
                outfile: output_file.to_string(),
                split_window: args.split_window,
            };

            merger.merge_windows();
        }
    }
}

// }
// }
