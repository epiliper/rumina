use bam::Record;
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

    // retrieve bundles from umi list
    let mut n = 0;
    let grouper: Grouper;

    for position in bottomhash.bottom_dict.keys() {
        for k in bottomhash.bottom_dict[position].keys() {
            let bundle = bottomhash.bottom_dict[position].get(k).unwrap();
            let umis = bundle
                .keys()
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let processor = processor::Processor { umis: &umis };

            let mut counts: HashMap<String, i32> = HashMap::new();
            for umi in &umis {
                // counts.insert(umi.to_string(), bundle[umi].count);
                counts.entry(umi.to_string()).or_insert(bundle[umi].count);
            }

            println! {"length of counts {:?}", counts.len()};

            let groupies = processor.main_grouper(counts.clone());
            n += 1;
            println! {"{}", n};


        }
    }

    let elapsed = now.elapsed();
    println!{"Time elapsed {:.2?}", elapsed};
}
