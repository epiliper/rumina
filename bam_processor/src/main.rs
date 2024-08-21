#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::fs::OpenOptions;
use crate::read_io::ChunkProcessor;
use crate::read_io::GroupReport;
use bam::bam_writer::BamWriterBuilder;
use bam::RecordWriter;
use clap::ValueEnum;
use indexmap::IndexMap;
use std::fs;
use std::fs::File;
use std::hash::DefaultHasher;
use std::hash::Hasher;
use std::io::Write;
use std::mem::drop;
use std::path::Path;
use std::time::Instant;

use bam::Record;

use std::hash::Hash;

use clap::Parser;
use rayon::ThreadPoolBuilder;

use parking_lot::Mutex;
use std::sync::Arc;

mod bottomhash;
mod deduplicator;
mod grouper;
mod read_io;
mod read_picker;

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
    #[arg(long = "length")]
    length: bool,
    #[arg(long = "only-group")]
    only_group: bool,
    #[arg(long = "singletons")]
    singletons: bool,
}

fn main() {
    let args = Args::parse();

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("ERROR: Invalid number of threads specified!");

    let now = Instant::now();

    // holds organized reads
    let bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };

    let input_file = args.input;
    let output_file = args.output;
    let separator = args.separator;
    let grouping_method = args.grouping_method;

    // create tag seed based on input file name
    let mut hasher = DefaultHasher::new();
    input_file.hash(&mut hasher);
    let seed = hasher.finish();

    let reads_to_write: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));

    let bam = bam::BamReader::from_path(&input_file, (args.threads) as u16).unwrap();
    let header = bam::BamReader::from_path(&input_file, 0)
        .unwrap()
        .header()
        .clone();

    let out_bam = output_file.clone();

    // records min and max reads per group
    let min_maxes: Arc<Mutex<GroupReport>> = Arc::new(Mutex::new(GroupReport {
        min_reads: i64::MAX,
        min_reads_group: *b"NONENONE",
        max_reads: 0,
        max_reads_group: *b"NONENONE",
        num_passing_groups: 0,
        num_groups: 0,
        num_umis: 0,
        num_reads_input_file: 0,
        num_reads_output_file: 0,
    }));

    // holds filtered reads awaiting writing to output bam file
    let mut read_handler = ChunkProcessor {
        separator: &separator,
        reads_to_output: Arc::clone(&reads_to_write),
        min_max: Arc::clone(&min_maxes),
        grouping_method,
        group_by_length: args.length,
        seed: seed,
        only_group: args.only_group,
        singletons: args.singletons,
        read_counter: 0,
    };

    // do grouping and processing
    read_handler.process_chunks(bam, bottomhash);
    let num_reads_in = read_handler.read_counter;

    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};

    let mut bam_writer = BamWriterBuilder::new()
        .additional_threads(args.threads.try_into().unwrap())
        .from_path(&out_bam, header)
        .unwrap();

    for read in reads_to_write.lock().drain(0..) {
        bam_writer.write(&read);
    }

    drop(read_handler);

    let mut group_report = Arc::try_unwrap(min_maxes).unwrap().into_inner();
    group_report.num_reads_input_file = num_reads_in;

    // report on min and max number of reads per group
    // this creates minmax.txt
    if group_report.min_reads != i64::MAX {
        println! {"minimum number of reads per group:     {},    group: {:?}", group_report.min_reads, String::from_utf8(group_report.min_reads_group.to_vec()).unwrap()};
        println! {"maximum number of reads per group:     {},    group: {:?}", group_report.max_reads, String::from_utf8(group_report.max_reads_group.to_vec()).unwrap()};
        println! {"writing remaining reads..."};

        let minmax_file = Path::new(&output_file).parent().unwrap().join("minmax.txt");

        if !minmax_file.exists() {
            File::create(&minmax_file);
        }
        let mut f = OpenOptions::new()
            .write(true)
            .append(true)
            .open(&minmax_file)
            .expect("unable to open minmax file");

        let _ = f.write(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                group_report.num_reads_input_file,
                group_report.num_reads_output_file,
                group_report.num_umis,
                group_report.num_groups,
                group_report.num_passing_groups,
                String::from_utf8(group_report.min_reads_group.to_vec()).unwrap(),
                group_report.min_reads,
                String::from_utf8(group_report.max_reads_group.to_vec()).unwrap(),
                group_report.max_reads,
            )
            .as_bytes(),
        );
    }
}
