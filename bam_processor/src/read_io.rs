use crate::bottomhash::BottomHashMap;
use crate::processor;
use crate::Grouper;
use bam::BamReader;
use bam::Record;
use indicatif::ProgressBar;
use std::collections::HashMap;
use std::fs::File;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}


pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub reads_to_output: &'a mut Vec<Record>,
    pub counter: i64,
    pub chunksize: i32,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(
        bottomhash: &mut BottomHashMap,
        mut reads_to_output: &mut Vec<Record>,
        grouper: & mut Grouper,
    ) {
        let max = bottomhash.bottom_dict.keys().len();
        let progressbar = ProgressBar::new(max.try_into().unwrap());

        for (i, position) in bottomhash.bottom_dict.values_mut().enumerate() {
            print! {"\rProcessing bundle {} of {}", i, max};
            let bundle = position.shift_remove(&0).unwrap();
            let umis = bundle
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
        progressbar.finish();
    }

pub fn pull_read(read: &Record, bottomhash: &mut BottomHashMap, separator: &String, mut counter: i64, chunksize: i32, mut reads_to_output: &mut Vec<Record>, grouper: & mut Grouper) { 
        if read.flag().is_mapped()
            && read.flag().is_reverse_strand()
        {
            if read.calculate_end() + 1 % chunksize == 0 {
            }
            bottomhash.update_dict(&(&read.calculate_end() + 1), 0, &get_umi(&read, separator), &read);
        // otherwise, use it's first position to reference
        } else if read.flag().is_mapped() {
            if read.start() + 1 % chunksize == 0 {
                Self::group_reads(bottomhash, reads_to_output, grouper);
            }
            bottomhash.update_dict(&(&read.start() + 1), 0, &get_umi(&read, separator), &read);
            }
    }

    pub fn process_chunks(& mut self, input_file: BamReader<File>, mut bottomhash: BottomHashMap, mut grouper: Grouper,) {

        let mut counter = 0;

        for r in input_file {
            let read = &r.unwrap();
            Self::pull_read(read, &mut bottomhash, self.separator, self.counter, self.chunksize, self.reads_to_output, &mut grouper);
            counter += 1;
            if counter % 100_000 == 0 {
                println! {"\rRead in {counter} reads" }
            }
        }

        Self::group_reads(&mut bottomhash, self.reads_to_output, & mut grouper);



    }
}

