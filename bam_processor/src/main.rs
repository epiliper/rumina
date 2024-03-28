use bam::Record;
use std::collections::HashMap;
use std::env;

mod bottomhash;

fn get_umi(record: &Record) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(":").last().unwrap().to_string()
}

fn main() {
    // let mut bottom_dict: PositionKey = HashMap::new();
    let mut master_dic = bottomhash::BottomHashMap {
        bottom_dict: HashMap::new(),
    };

    let args: Vec<String> = env::args().collect();
    let input_file = &args[1];
    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();
    let mut n: i64 = 0;

    // implement get_bundles()
    for read in bam {
        if read.as_ref().unwrap().flag().is_reverse_strand()
            | !read.as_ref().unwrap().flag().is_mapped()
        // if read.unwrap().flag().is_reverse_strand()
        {
            continue;
        } else {
            let r1 = read.as_ref().unwrap();

            master_dic.update_dict(&r1.start(), 0, &get_umi(&r1), &r1);
            n += 1;

            //print the millionth iteration dict for debugging
            if n == 1_000_000 {
                println!("{:?}", &master_dic.bottom_dict);
            }

            if n % 100_000 == 0 {
                println! {"Processed {n} reads" }
            }
        }
    }
}
