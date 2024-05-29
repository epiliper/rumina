use crate::grouper::Grouper;
use std::path::Path;
use crate::read_io::ChunkProcessor;
use bam::bam_writer::BamWriterBuilder;
use bam::Record;
use bam::RecordWriter;
use indexmap::IndexMap;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelDrainRange; use rayon::iter::ParallelIterator;
use std::env;
use std::fs;
use std::io;
use std::io::Write;
use std::time::Instant;
use std::fs::File;
use polars::prelude::*;
use std::mem;
use crate::fs::OpenOptions;

use parking_lot::Mutex;
use std::sync::Arc;

mod bottomhash;
mod grouper;
mod processor;
mod read_io;

fn main() {
    let args: Vec<String> = env::args().collect();
    let report_list: Arc<Mutex<Vec<DataFrame>>> = Arc::new(Mutex::new(Vec::new()));

   let now = Instant::now();

    let bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };
    let input_file = &args[1];
    let output_file = &args[2];
    let separator = &args[3];

    // get info on min and max number of reads per group
   
    println! {"using separator = {}", separator};
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

    let min_maxes: Arc<Mutex<Vec<i64>>> = Arc::new(Mutex::new(Vec::new()));
    let reads_to_spit: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));

    let mut read_handler = ChunkProcessor {
        separator: separator,
        reads_to_output: Arc::clone(&reads_to_spit),
        reports: Arc::clone(&report_list),
        min_max: Arc::clone(&min_maxes),
    };

    read_handler.process_chunks(bam, bottomhash);

    println! {"Writing {} reads...", reads_to_spit.lock().len()};
    for read in reads_to_spit.lock().iter() {
        outfile.write(read).unwrap();
    }

    println! {"Writing results to .csv..."};

    mem::drop(read_handler);

    // let mut report = Arc::into_inner(report_list).unwrap().into_inner();
    // let chunksize = report.len() / 8;

    // let mut new_report = DataFrame::default();

    // let subframes: Arc<Mutex<Vec<DataFrame>>> = Arc::new(Mutex::new(Vec::new()));
    // report.par_drain(0..).chunks(chunksize).for_each(
    //     |x| {
    //         let mut miniframe = DataFrame::default();
    //         for frame in x {
    //             miniframe.vstack_mut(&frame).unwrap();
    //         }
    //         subframes.lock().push(miniframe);
    //     }
    // );

    // subframes.lock().iter().for_each(
    //     |x| {
    //         new_report.vstack_mut(x).unwrap();
    //     }
    // );

    // new_report.align_chunks();

    // let report_name = format!("{}{}", input_file.split(".bam").next().unwrap(), "_report.csv");
    // let mut file = File::create(report_name).expect("Could not create file!");
    // CsvWriter::new(&mut file)
    //     .n_threads(8)
    //     .finish(&mut new_report)
    //     .unwrap();

    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};

    let min_maxes = Arc::into_inner(min_maxes).unwrap().into_inner();


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




        // fs::write(&Path::new(input_file).parent().unwrap().join("minmax.txt"), format!("{}\t{}", min, max).as_bytes()).unwrap();

    }




}
