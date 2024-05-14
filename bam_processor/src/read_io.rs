use crate::bottomhash::BottomHashMap;
use crate::processor;
use crate::Grouper;
use bam::BamReader;
use bam::Record;
use indicatif::ProgressBar;
use parking_lot::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::sync::Arc;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub reads_to_output: Arc<Mutex<Vec<Record>>>,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(&mut self, bottomhash: &mut BottomHashMap) {
        let progressbar = ProgressBar::new(bottomhash.bottom_dict.keys().len().try_into().unwrap());
        bottomhash.bottom_dict.par_values_mut().for_each(|x| {
            let bundle = x.shift_remove(&0).unwrap();
            let umis = bundle
                .keys()
                .map(|x| x.to_string())
                .collect::<Vec<String>>();

            let processor = processor::Processor { umis: &umis };

            let mut counts: HashMap<&String, i32> = HashMap::new();

            for umi in &umis {
                counts.entry(umi).or_insert(bundle[umi].count);
            }

            let mut grouper = Grouper {};
            let groupies = processor.main_grouper(counts);
            let tagged_reads = grouper.tag_records(groupies, bundle);
            let mut out = self.reads_to_output.lock();

            match tagged_reads {
                Some(tagged_reads) => {
                    out.extend(tagged_reads);
                    drop(out);
                }
                None => (),
            }
            progressbar.inc(1);
        });
        progressbar.finish();
    }

    pub fn pull_read(&mut self, read: &Record, bottomhash: &mut BottomHashMap, separator: &String) {
        // if read is reverse to reference, group it by its last aligned base to the reference
        if read.flag().is_mapped() && read.flag().is_reverse_strand() {
            bottomhash.update_dict(
                &(&read.calculate_end() + 1),
                0,
                &get_umi(&read, separator),
                &read,
            );

        // otherwise, use its first position to reference
        } else if read.flag().is_mapped() {
            bottomhash.update_dict(&(&read.start() + 1), 0, &get_umi(&read, separator), &read);
        }
    }

    pub fn process_chunks(&mut self, input_file: BamReader<File>, mut bottomhash: BottomHashMap) {
        let mut counter = 0;

        for r in input_file {
            let read = &r.unwrap();
            Self::pull_read(self, read, &mut bottomhash, self.separator);

            counter += 1;
            if counter % 100_000 == 0 {
                print! {"\rRead in {counter} reads" }
            }
        }
        print! {"\r Grouping {counter} reads...\n"}

        Self::group_reads(self, &mut bottomhash);
    }
}
