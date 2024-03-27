extern crate bam;
extern crate rust_htslib;

use bam::Record;
use std::collections::HashMap;
use std::env;

fn get_umi(record: &Record) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(":").last().unwrap().to_string()
}

// the data structure used to cluster UMIs. These are all inserted into a bottom_level dictionary.
// top key: PositionKey
// bottom key: ReadsCount
type PositionKey = HashMap<i32, KeyUMI>; //every position has a key
type KeyUMI = HashMap<i32, UMIReads>; // every key has a UMI
type UMIReads = HashMap<String, ReadsAndCount>; // every UMI has a set of reads

#[derive(Debug)]
pub struct ReadsAndCount {
    pub reads: Vec<Record>,
    pub count: i32
}

impl ReadsAndCount {
    fn up(&mut self) {
        self.count += 1;
    }
}

fn update_dict(
    position: &i32,
    key: i32,
    umi: &String,
    read: &Record,
    // read_count: i32,
    final_dict: &mut PositionKey,
) {
    // at every position, at each key, for each umi, append a bam record to list of reads
    final_dict
        .entry((*position).into())
        .or_default()
        .entry((key).into())
        .or_default()
        .entry(umi.into())
        .or_insert(ReadsAndCount {
            reads: Vec::new(),
            count: 0
        })
            .reads.push(read.clone());

    // for every record appended, increment the count field by 1
    final_dict
        .entry((*position).into())
        .or_default()
        .entry((key).into())
        .or_default()
        .entry(umi.into())
        .or_insert(ReadsAndCount {
            reads: Vec::new(),
            count: 0
        })
            .up();
}

fn main() {
    let mut bottom_dict: PositionKey = HashMap::new();

    let args: Vec<String> = env::args().collect();
    let input_file = &args[1];
    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();
    let mut n: i64 = 0;

    for read in bam {
        if read.as_ref().unwrap().flag().is_reverse_strand()
        // if read.unwrap().flag().is_reverse_strand()
        {
            continue;
        } else {
            let r1 = read.as_ref().unwrap();

            update_dict(&r1.start(), 0, &get_umi(&r1), &r1, &mut bottom_dict);
            n += 1;
            
            //print the millionth iteration dict for debugging
            if n == 1_000_000 {
                println!("{:?}", &bottom_dict);
            }

            if n % 100_000 == 0 {
                println!{"Processed {n} reads" }
            } 
        }

    }
}
