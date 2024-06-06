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
use crate::GroupingMethod;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

#[derive(Default, Debug)]
pub struct GroupReport {
    pub min_reads: i64,
    pub max_reads: i64,
    pub min_reads_group: [u8; 8],
    pub max_reads_group: [u8; 8],
    pub num_passing_groups: i64,
    pub num_groups: i64,
}

pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub reads_to_output: Arc<Mutex<Vec<Record>>>,
    pub min_max: Arc<Mutex<GroupReport>>,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(&mut self, bottomhash: &mut BottomHashMap, grouping_method: GroupingMethod) {
        let progressbar = ProgressBar::new(bottomhash.bottom_dict.keys().len().try_into().unwrap());
        bottomhash.bottom_dict.par_values_mut().for_each(|x| {
            let bundle = x.shift_remove(&0).unwrap();
            let umis = bundle
                .keys()
                .map(|x| x.to_string())
                .collect::<Vec<String>>();

            let processor = processor::Processor { umis: &umis };
            let mut counts: HashMap<&String, i32> = HashMap::new();

            // get number of reads for each raw UMI
            for umi in &umis {
                counts.entry(umi).or_insert(bundle[umi].count);
            }

            let groupies: (HashMap<&String, i32>, Option<Vec<Vec<&String>>>);
            let mut grouper = Grouper {};

            match grouping_method {
                GroupingMethod::Directional => {groupies = processor.directional_clustering(counts);},
                GroupingMethod::Raw => {groupies = processor.no_clustering(counts);}
             }
            let tagged_reads = grouper.tag_records(groupies, bundle);

            let mut out = self.reads_to_output.lock();
            let mut min_max = self.min_max.lock();

            // update the groups with mininum and maximum observed reads
            match tagged_reads {
                Some(tagged_reads) => {
                    out.extend(tagged_reads.1);
                    drop(out);

                    match tagged_reads.0 {
                        Some(x) => {

                            if x.max_reads > min_max.max_reads {
                                min_max.max_reads = x.max_reads;
                                min_max.max_reads_group = x.max_reads_group;
                            }

                            if x.min_reads < min_max.min_reads {
                                min_max.min_reads = x.min_reads;
                                min_max.min_reads_group = x.min_reads_group;
                            }

                            // count the number of UMI groups used in consensus
                            min_max.num_passing_groups += x.num_passing_groups;
                            min_max.num_groups += x.num_groups;

                        }
                        _ => ()
                    }
                }
                None => (),
            }
            progressbar.inc(1);
        });
        progressbar.finish();
    }

    // organize reads in bottomhash based on position
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

    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(&mut self, input_file: BamReader<File>, mut bottomhash: BottomHashMap, grouping_method: GroupingMethod) {
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

        Self::group_reads(self, &mut bottomhash, grouping_method);
    }

}
