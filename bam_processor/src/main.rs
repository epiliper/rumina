use crate::fs::OpenOptions;
use crate::read_io::ChunkProcessor;
use crate::read_io::GroupReport;
use bam::BamWriter;
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
use std::sync::mpsc::channel;
use std::thread;

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
    Bidirectional,
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

    // create tag seed based on input file name
    let mut hasher = DefaultHasher::new();
    input_file.hash(&mut hasher);
    let seed = hasher.finish();

    let bam = bam::BamReader::from_path(&input_file, (args.threads - 1) as u16).unwrap();
    let header = bam::BamReader::from_path(&input_file, 0)
        .unwrap()
        .header()
        .clone();

    // This gets the reads processed per position and sends them for file writing
    // this is pending renaming/polish
    let (tx, rx) = channel::<Vec<Record>>();

    let out_bam = output_file.clone();

    let _ = BamWriter::from_path(&out_bam, header.clone()).unwrap();

    // writes reads in the read channel to output bam
    let writer_handle = thread::spawn(move || {
        let mut bam_writer = BamWriter::from_path(out_bam, header).unwrap();

        while let Ok(reads) = rx.recv() {
            for read in reads {
                bam_writer.write(&read).unwrap();
            }
            bam_writer.flush().unwrap();
        }
    });

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
    let mut read_handler = ChunkProcessor {
        separator: &separator,
        reads_to_output: tx,
        min_max: Arc::clone(&min_maxes),
        grouping_method,
        group_by_length: args.length,
        seed: seed,
    };

    // do grouping and processing
    read_handler.process_chunks(bam, bottomhash);

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

    writer_handle.join().unwrap();
}
