extern crate bam;
extern crate rust_htslib;

use bam::record::tags::TagValue;
use std::collections::HashMap;

use std::fs::OpenOptions;
use std::io::prelude::*;

use std::env;

fn main() {

    let args: Vec<String> = env::args().collect();

    let input_file = &args[1];
    println!("{}", input_file);

    let bam = bam::BamReader::from_path(&input_file, 4).unwrap();

    let mut umis = HashMap::new();

    let threshold: i64 = 1;

    println!("Identifying onesies...");

    for read in bam {
        // let record = read.unwrap();
        match read.unwrap().tags().get(b"UG") {
            Some(TagValue::Int(value, _u32)) => {
                umis.entry(value).and_modify(|counter| * counter += 1).or_insert(1);
            }
            _ => panic!("Unable to retrieve UG tags. Check that all reads have UG tags.")
        }
    };

    let mut _blacklist_file_name = input_file.split(".bam").next().unwrap().to_string() + "_onesies.txt";

    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .append(true)
        .open(_blacklist_file_name)
        .unwrap();

    for (k, _v) in umis.iter()
        .filter(|&(_k, v)| v > &threshold) {
        writeln!(file, "{}", k).expect("Error writing blacklist file");
    }

    println!("Onesies identified! Purge commencing...")

}
