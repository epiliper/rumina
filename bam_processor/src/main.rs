use bam::Record;
use std::collections::HashMap;
use std::env;

mod bottomhash;
mod processor;

fn get_umi(record: &Record) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(":").last().unwrap().to_string()
}

fn main() {
    let mut bottomhash = bottomhash::BottomHashMap {
        bottom_dict: HashMap::new(),
    };
    let args: Vec<String> = env::args().collect();
    let input_file = &args[1];
    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();
    let mut n: i64 = 0;

    for read in bam {
        // don't process unmappped or reverse strands.
        if read.as_ref().unwrap().flag().is_reverse_strand()
            | !read.as_ref().unwrap().flag().is_mapped()
        {
            continue;
        } else {
            let r1 = read.as_ref().unwrap();

            // get the Record itself, plus its UMI
            bottomhash.update_dict(&r1.start(), 0, &get_umi(&r1), &r1);
            n += 1;

            if n % 100_000 == 0 {
                println! {"Read in {n} reads" }
            }
        }
    }

    let mut counts: HashMap<String, i32> = HashMap::new();

    let processor = processor::Processor{};

    // retrieve bundles from umi list
    for (bundle, key) in bottomhash.iter() {
        let umis = bundle.keys().map(|x| x.to_string()).collect::<Vec<String>>();

        for umi in &umis {
            // counts.insert(umi.to_string(), bundle[umi].count);
            counts.entry(umi.to_string()).or_insert(bundle[umi].count);
        }
        let groupies = processor.main_grouper(umis, counts.clone());
        println!{"{:?}", groupies};
    }




}
