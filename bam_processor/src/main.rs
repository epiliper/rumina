use crate::fs::OpenOptions;
use crate::read_io::ChunkProcessor;
use crate::read_io::GroupReport;
use bam::bam_writer::BamWriterBuilder;
use bam::Record;
use bam::RecordWriter;
use clap::ValueEnum;
use indexmap::IndexMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::mem::drop;
use std::path::Path;
use std::time::Instant;

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
}

fn main() {
    let args = Args::parse();

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
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

    let bam = bam::BamReader::from_path(&input_file, (args.threads - 1) as u16).unwrap();
    let header = bam::BamReader::from_path(&input_file, 0)
        .unwrap()
        .header()
        .clone();

    let mut outfile = BamWriterBuilder::from_path(
        &mut BamWriterBuilder::new().additional_threads((args.threads - 1) as u16),
        &output_file,
        header,
    )
    .unwrap();

    // records min and max reads per group
    let min_maxes: Arc<Mutex<GroupReport>> = Arc::new(Mutex::new(GroupReport {
        min_reads: i64::MAX,
        min_reads_group: *b"NONENONE",
        max_reads: 0,
        max_reads_group: *b"NONENONE",
        num_passing_groups: 0,
        num_groups: 0,
        num_umis: 0,
    }));

    // holds filtered reads awaiting writing to output bam file
    let reads_to_spit: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));

    let mut read_handler = ChunkProcessor {
        separator: &separator,
        reads_to_output: Arc::clone(&reads_to_spit),
        min_max: Arc::clone(&min_maxes),
        grouping_method,
    };

    // do grouping and processing
    read_handler.process_chunks(bam, bottomhash);

    // write final reads to output
    println! {"Writing {} reads...", reads_to_spit.lock().len()};
    for read in reads_to_spit.lock().iter() {
        outfile.write(read).unwrap();
    }

    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};

    drop(read_handler);
    let group_report = Arc::try_unwrap(min_maxes).unwrap().into_inner();

    // report on min and max number of reads per group
    // this creates minmax.txt
    if group_report.min_reads != i64::MAX {
        println! {"minimum number of reads per group:     {},    group: {:?}", group_report.min_reads, String::from_utf8(group_report.min_reads_group.to_vec()).unwrap()};
        println! {"maximum number of reads per group:     {},    group: {:?}", group_report.max_reads, String::from_utf8(group_report.max_reads_group.to_vec()).unwrap()};

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
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                String::from_utf8(group_report.min_reads_group.to_vec()).unwrap(),
                group_report.min_reads,
                String::from_utf8(group_report.max_reads_group.to_vec()).unwrap(),
                group_report.max_reads,
                group_report.num_passing_groups,
                group_report.num_groups,
                group_report.num_umis,
            )
            .as_bytes(),
        );
    }
}
