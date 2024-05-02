use crate::grouper::Grouper;
use crate::read_io::ChunkProcessor;
use bam::BamWriter;
use bam::Record;
use bam::RecordWriter;
use indexmap::IndexMap;
use std::env;
use std::time::Instant;

mod bottomhash;
mod grouper;
mod processor;
mod read_io;

fn main() {
    let now = Instant::now();

    let mut bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };
    let args: Vec<String> = env::args().collect();
    println! {"using separator = {}", args[3]};
    let input_file = &args[1];
    let output_file = &args[2];
    let separator = &args[3];
    let bam = bam::BamReader::from_path(&input_file, 6).unwrap();
    let header = bam::BamReader::from_path(&input_file, 6)
        .unwrap()
        .header()
        .clone();

    let n: i64 = 0;
    let mut outfile = BamWriter::from_path(&output_file, header).unwrap();
    let grouper = Grouper { num: 0 };
    let mut reads_to_spit: Vec<Record> = Vec::new();


    let mut read_handler = ChunkProcessor{
        separator: separator, 
        reads_to_output: &mut reads_to_spit, 
        counter: n, 
        chunksize: 500
    };

    read_handler.process_chunks(bam, bottomhash, grouper);

    println! {"Writing {} reads to {}", &reads_to_spit.len(), output_file};

    reads_to_spit.iter().for_each(|x| outfile.write(x).unwrap());

    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};
}
