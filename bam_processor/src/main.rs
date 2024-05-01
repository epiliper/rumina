use crate::grouper::Grouper;
use crate::read_io::group_reads;
use bam::BamWriter;
use bam::Record;
use bam::RecordWriter;
use indexmap::IndexMap;
use std::env;
use std::time::Instant;
use crate::read_io::pull_reads;

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
    let bam = bam::BamReader::from_path(&input_file, 6).unwrap();
    let header = bam::BamReader::from_path(&input_file, 6)
        .unwrap()
        .header()
        .clone();

    let n: i64 = 0;
    let mut outfile = BamWriter::from_path(&output_file, header).unwrap();
    let grouper = Grouper { num: 0 };
    let mut reads_to_spit: Vec<Record> = Vec::new();

    pull_reads(bam, & mut bottomhash, &args[3], n);
    group_reads(& mut bottomhash, & mut reads_to_spit, grouper);


    // progressbar.finish();
    println! {"Writing {} reads to {}", &reads_to_spit.len(), output_file};

    reads_to_spit.iter().for_each(|x| outfile.write(x).unwrap());

    let elapsed = now.elapsed();
    println! {"Time elapsed {:.2?}", elapsed};
}
