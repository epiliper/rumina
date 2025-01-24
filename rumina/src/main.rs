#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::args::{parse_args, GroupingMethod};
use crate::group_report::GroupReport;
use crate::main_dedup::{init_processor, process_chunks};
use crate::pair_merger::PairMerger;
use crate::utils::index_bam;
use colored::Colorize;
use indexmap::IndexMap;
use std::fs::{create_dir, read_dir, remove_file};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use utils::{gen_outfile_name, get_file_ext};

use log::{error, info, LevelFilter};
use rayon::ThreadPoolBuilder;

use parking_lot::Mutex;
use std::collections::HashMap;
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

    for input_file in infiles {
        // create tag seed based on input file name
        let output_file = gen_outfile_name(Some(&args.outdir), "RUMINA", &input_file.1);

        let mut hasher = DefaultHasher::new();
        input_file.hash(&mut hasher);
        let seed = hasher.finish();

        // create deduplication report
        let min_maxes: Arc<Mutex<GroupReport>> = Arc::new(Mutex::new(GroupReport::new()));

        let (bam_reader, bam_writer, mut read_handler) = init_processor(
            input_file.0.clone(),
            output_file.to_string(),
            args.grouping_method.clone(),
            args.threads,
            args.split_window,
            args.length,
            args.only_group,
            args.singletons,
            args.r1_only,
            min_maxes.clone(),
            seed,
        );

        info!("{:?}", read_handler);

        // holds filtered reads awaiting writing to output bam file
        // do grouping and processing
        process_chunks(&mut read_handler, bam_reader, &args.separator, bam_writer);
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

        index_bam(&output_file, args.threads).unwrap();

        if let Some(ref ref_fasta) = args.merge_pairs {
            let mut merger = PairMerger {
                ref_fasta: ref_fasta.to_string(),
                min_overlap_bp: args.min_overlap_bp,
                threads: args.threads,
                infile: output_file.to_string(),
                outfile: gen_outfile_name(None, "MERGED", &output_file),
                split_window: args.split_window,
            };

            info!("{:?}", merger);

            let merge_report = merger.merge_windows();
            remove_file(output_file).unwrap();
            print!("{merge_report}");
        }
    }
}
