use bam::Record;
use std::collections::HashSet;
use bam::BamWriter;
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
    let bam = bam::BamReader::from_path(&input_file, 6).unwrap();
    let header = bam::BamReader::from_path(&input_file, 6).unwrap().header().clone();
    let mut n: i64 = 0;


    for read in bam {
        // don't process unmappped or reverse strands.
        if read.as_ref().unwrap().flag().is_reverse_strand()
            | !(read.as_ref().unwrap().flag().is_mapped())
        {
            continue;
        } else {
            let r1 = read.as_ref().unwrap();
            if get_umi(r1).len() != 11 {
                continue;
            } else {
                // get the Record itself, plus its UMI
                bottomhash.update_dict(&r1.start(), 0, &get_umi(&r1), &r1);
    
                n += 1;

                if n % 100_000 == 0 {
                    println! {"Read in {n} reads" }
                }
            }
        }
    }

    let outfile = BamWriter::from_path("test/output.bam", header).unwrap();

    let grouper = Grouper{};

    let mut reads_to_spit: Vec<Record> = Vec::new();


    for position in bottomhash.bottom_dict.values_mut() {
        let bundle = position.shift_remove(&0).unwrap();
        let mut umis = bundle.keys().map(|x| x.to_string()).collect::<Vec<String>>();
        umis.dedup();

        let processor = processor::Processor {umis: &umis};

        let mut counts: HashMap<&String, i32> = HashMap::new();

        for umi in &umis {
            counts.entry(umi).or_insert(bundle[umi].count);

        }

        println!{"length of counts {:?}", counts.len()};

        let groupies = processor.main_grouper(counts);

        grouper.tag_records(groupies, bundle, & mut reads_to_spit);



    }
    

    

        
    let elapsed = now.elapsed();
    println!{"Time elapsed {:.2?}", elapsed};
}
