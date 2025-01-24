#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, Command, GroupingMethod};
use crate::group_report::GroupReport;
use crate::main_dedup::{init_processor, process_chunks};
use crate::pair_merger::PairMerger;
use crate::utils::index_bam;
use colored::Colorize;
use indexmap::IndexMap;
use std::hash::{DefaultHasher, Hash, Hasher};

use log::{info, LevelFilter};
use rayon::ThreadPoolBuilder;

use parking_lot::Mutex;
use std::sync::Arc;

mod args;
mod bottomhash;
mod deduplicator;
mod group_report;
mod grouper;
mod main_dedup;
mod merge;
mod merge_report;
mod pair_merger;
mod progbars;
mod read_picker;
mod readkey;
mod realign;
mod utils;
mod window_processor;

fn main() {
    let args = parse_args();
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

            simple_logging::log_to_file("rumina_group.log", LevelFilter::Info)
                .expect("Failed to create log for grouping command!");

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

            info!("{:?}", read_handler);

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
            if !group_report.is_blank() {
                println!("{}", "DONE".green());

                group_report.write_to_report_file(&output_file);
                println!("{}", group_report);
            }

            index_bam(output_file, args.threads).unwrap();
        }

        Command::Merge {
            ref_fasta,
            min_overlap_bp,
        } => {
            simple_logging::log_to_file("rumina_merge.log", LevelFilter::Info)
                .expect("Failed to create log for merge command!");

            let mut merger = PairMerger {
                ref_fasta,
                min_overlap_bp,
                threads: args.threads,
                infile: input_file.to_string(),
                outfile: output_file.to_string(),
                split_window: args.split_window,
            };

            info!("{:?}", merger);

            let merge_report = merger.merge_windows();
            print!("{merge_report}");
        }
    }
}
