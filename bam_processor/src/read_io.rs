use crate::bottomhash::BottomHashMap;
use indicatif::ProgressBar;
use crate::Grouper;
use crate::processor;
use std::collections::HashMap;
use bam::Record;
use bam::BamReader;
use std::fs::File;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

pub fn pull_reads(bam:BamReader<File>, bottomhash: & mut BottomHashMap, separator:&String, mut counter: i64) {

    for read in bam {
        if read.as_ref().unwrap().flag().is_mapped() && read.as_ref().unwrap().flag().is_reverse_strand(){
            let r1 = &read.unwrap();
            bottomhash.update_dict(&(&r1.calculate_end() + 1), 0, &get_umi(&r1, separator), &r1);
            counter += 1;
            if counter % 100_000 == 0 {
                println! {"\rRead in {counter} reads" }
            }
        } else if read.as_ref().unwrap().flag().is_mapped() {
            let r1 = &read.unwrap();
            bottomhash.update_dict(&(&r1.start() + 1), 0, &get_umi(&r1, separator), &r1);
            counter += 1;
            if counter % 100_000 == 0 {
                println! {"\rRead in {counter} reads" }
            }
        }
    }
}

pub fn group_reads(bottomhash: & mut BottomHashMap, mut reads_to_output: & mut Vec<Record>, mut grouper: Grouper) {
        let max = bottomhash.bottom_dict.keys().len();
        let progressbar = ProgressBar::new(max.try_into().unwrap());

        for (i, position) in bottomhash.bottom_dict.values_mut().enumerate() {
        print! {"\rProcessing bundle {} of {}", i, max};
        let bundle = position.shift_remove(&0).unwrap();
        let mut umis = bundle
            .keys()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let processor = processor::Processor { umis: &umis };

        let mut counts: HashMap<&String, i32> = HashMap::new();

        for umi in &umis {
            counts.entry(umi).or_insert(bundle[umi].count);
        }

        let groupies = processor.main_grouper(counts);
        grouper.tag_records(groupies, bundle, &mut reads_to_output);
        progressbar.inc(1);

    }


}

