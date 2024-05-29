use crate::grouper::Grouper;
use std::path::Path;
use crate::read_io::ChunkProcessor;
use bam::bam_writer::BamWriterBuilder;
use bam::Record;
use bam::RecordWriter;
use indexmap::IndexMap;
use std::env;
use std::fs;
use std::io::Write;
use std::time::Instant;
use std::fs::File;
use std::mem::drop;
use crate::fs::OpenOptions;

use parking_lot::Mutex;
use std::sync::Arc;

mod bottomhash;
mod grouper;
mod processor;
mod read_io;
mod dedup_correct;

fn main() {
    let args: Vec<String> = env::args().collect();

   let now = Instant::now();

    // holds organized reads
    let bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };
    let input_file = &args[1];
    let output_file = &args[2];
    let separator = &args[3];

    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();
    let header = bam::BamReader::from_path(&input_file, 0)
        .unwrap()
        .header()
        .clone();

    let mut outfile = BamWriterBuilder::from_path(
        &mut BamWriterBuilder::new().additional_threads(5),
        &output_file,
        header,
    )
    .unwrap();

    // records min and max reads per group
    let min_maxes: Arc<Mutex<Vec<i64>>> = Arc::new(Mutex::new(Vec::new()));

    // holds filtered reads awaiting writing to output bam file
    let reads_to_spit: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));

    let mut read_handler = ChunkProcessor {
        separator: separator,
        reads_to_output: Arc::clone(&reads_to_spit),
        min_max: Arc::clone(&min_maxes),
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
    let min_maxes = Arc::try_unwrap(min_maxes).unwrap().into_inner();

    // report on min and max number of reads per group
    // this creates minmax.txt
    if !min_maxes.is_empty() {

        let min = &min_maxes.iter().min().unwrap();
        let max = &min_maxes.iter().max().unwrap();

        println!{"minimum number of reads per group     {}", min};
        println!{"maxmimum number of reads per group     {}", max};

        let minmax_file = Path::new(&output_file).parent().unwrap().join("minmax.txt");

        if !minmax_file.exists() {
            File::create(&minmax_file);
        }
        let mut f = OpenOptions::new()
            .write(true)
            .append(true)
            .open(&minmax_file)
            .expect("unable to open minmax file");

        f.write(format!("{}\t{}\n", min, max).as_bytes());

    }




}
