use bam::Record;
use bam::BamWriter;
use bam::RecordWriter;
use indexmap::IndexMap;
use std::collections::HashMap;
use std::env;
use std::time::Instant;
use crate::grouper::Grouper;

mod bottomhash;
mod processor;
mod grouper;

fn get_umi(record: &Record) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(":").last().unwrap().to_string()
}

fn main() { let now = Instant::now();

    let mut bottomhash = bottomhash::BottomHashMap {
        bottom_dict: IndexMap::new(),
    };
    let args: Vec<String> = env::args().collect();
    let input_file = &args[1];
    let output_file = &args[2];
    let bam = bam::BamReader::from_path(&input_file, 6).unwrap();
    let header = bam::BamReader::from_path(&input_file, 6).unwrap().header().clone();
    let mut n: i64 = 0;

    for read in bam {
        if read.as_ref().unwrap().flag().is_paired() {

            if read.as_ref().unwrap().flag().is_mapped() && read.as_ref().unwrap().flag().first_in_pair() && get_umi(read.as_ref().unwrap()).len() == 11 {
                let r1 = &read.unwrap();
                bottomhash.update_dict(&r1.start(), 0, &get_umi(&r1), &r1);
                n += 1;

                if n % 100_000 == 0 {
                    println!{"Read in {n} reads" }
                }

            }
        } else {
            if read.as_ref().unwrap().flag().is_mapped() && get_umi(read.as_ref().unwrap()).len() == 11 {
                let r1 = &read.unwrap();
                bottomhash.update_dict(&r1.start(), 0, &get_umi(&r1), &r1);
                n += 1;

                if n % 100_000 == 0 {
                    println!{"Read in {n} reads" }
                }

            }
        }
    }


    let mut outfile = BamWriter::from_path(&output_file, header).unwrap();

    let grouper = Grouper{};
    let max = bottomhash.bottom_dict.keys().len();

    let mut reads_to_spit: Vec<Record> = Vec::new();

    let mut counter = 1;

    for position in bottomhash.bottom_dict.values_mut() {
        print!{"\rProcessing bundle {} of {}", counter, max};
        let bundle = position.shift_remove(&0).unwrap();
        let mut umis = bundle.keys().map(|x| x.to_string()).collect::<Vec<String>>();
        umis.dedup();


        let processor = processor::Processor {umis: &umis};

        let mut counts: HashMap<&String, i32> = HashMap::new();

        for umi in &umis {
            counts.entry(umi).or_insert(bundle[umi].count);
        }


        let groupies = processor.main_grouper(counts);
        grouper.tag_records(groupies, bundle, & mut reads_to_spit);
        counter += 1;
    }

    println!{"Writing {} reads to {}", reads_to_spit.len(), output_file};
    
    reads_to_spit.iter().for_each(
        |x| outfile.write(x).unwrap()
    );

    let elapsed = now.elapsed();
    println!{"Time elapsed {:.2?}", elapsed};
}
