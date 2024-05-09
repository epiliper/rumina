use crate::grouper::Grouper;
use crate::read_io::ChunkProcessor;
use bam::bam_writer::BamWriterBuilder;
use bam::RecordWriter;
use bam::Record;
use indexmap::IndexMap;
use std::env;
use std::time::Instant;

use std::sync::{Arc};
use parking_lot::Mutex;

mod bottomhash;
mod grouper;
mod processor;
mod read_io;


fn main() {
    let args: Vec<String> = env::args().collect();

    // let pool = rayon::ThreadPoolBuilder::new();
    let now = Instant::now();
    // pool.num_threads(4);


    let bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };
    println! {"using separator = {}", args[3]};
    let input_file = &args[1];
    let output_file = &args[2];
    let separator = &args[3];

    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();
    let header = bam::BamReader::from_path(&input_file, 0)
        .unwrap()
        .header()
        .clone();

    let mut outfile = BamWriterBuilder::from_path(&mut BamWriterBuilder::new().additional_threads(5), &output_file, header).unwrap();
    let reads_to_spit: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));



    let mut read_handler = ChunkProcessor{
        separator: separator, 
        reads_to_output: Arc::clone(&reads_to_spit),
    };

    read_handler.process_chunks(bam, bottomhash);

    println!{"Writing {} reads...", reads_to_spit.lock().len()};
    for read in reads_to_spit.lock().iter() {
        outfile.write(read).unwrap();
    }


    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};
}
