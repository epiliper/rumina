use crate::bottomhash::BottomHashMap;
use crate::deduplicator::Deduplicator;
use crate::grouper::Grouper;
use crate::GroupingMethod;
use bam::BamReader;
use bam::Record;
use indicatif::ProgressBar;
use parking_lot::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::sync::mpsc::Sender;
use std::sync::Arc;

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
    pub num_umis: i64,
}

pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub reads_to_output: Sender<Vec<Record>>,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub only_group: bool,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam

    pub fn group_reads(&mut self, bottomhash: &mut BottomHashMap) {
        let progressbar = ProgressBar::new(bottomhash.bottom_dict.keys().len().try_into().unwrap());

        bottomhash.bottom_dict.par_drain(0..).for_each(|position| {
            for umi in position.1 {
                let mut umis_reads = umi.1;

                // sort UMIs by read count
                // note that this is an unstable sort, so we need to identify read-tied groups
                umis_reads.par_sort_unstable_by(|_umi1, count1, _umi2, count2| {
                    count2.count.cmp(&count1.count)
                });

                let umis = umis_reads
                    .keys()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();

                let processor = Grouper { umis: &umis };
                let mut counts: HashMap<&String, i32> = HashMap::with_capacity(umis_reads.len());

                // get number of reads for each raw UMI
                let mut num_umis = 0;
                for umi in &umis {
                    counts.entry(umi).or_insert(umis_reads[umi].count);
                    num_umis += 1;
                }

                let groupies: (HashMap<&String, i32>, Option<Vec<Vec<&String>>>);
                let mut grouper = Deduplicator {
                    seed: self.seed + position.0 as u64, // make seed unique per position
                    group_only: self.only_group,
                };

                // get grouping method
                match self.grouping_method {
                    GroupingMethod::Directional => {
                        groupies = processor.directional_clustering(counts);
                    }
                    GroupingMethod::Raw => {
                        groupies = processor.no_clustering(counts);
                    }
                    GroupingMethod::Acyclic => {
                        groupies = processor.bidirectional_clustering(counts);
                    }
                }
                let tagged_reads = grouper.tag_records(groupies, umis_reads);
                let mut min_max = self.min_max.lock();

                // update the groups with mininum and maximum observed reads
                match tagged_reads {
                    Some(tagged_reads) => {
                        self.reads_to_output
                            .send(tagged_reads.1)
                            .expect("Read channel not recieving!");

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
                                min_max.num_umis += num_umis;
                            }
                            _ => (),
                        }
                    }
                    None => (),
                }
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
    pub fn pull_read_w_length(
        &mut self,
        read: &Record,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        if read.flag().is_mapped() && read.flag().is_reverse_strand() {
            bottomhash.update_dict(
                &(&read.calculate_end() + 1),
                read.query_len() as i32,
                &get_umi(&read, separator),
                &read,
            );

        // otherwise, use its first position to reference
        } else if read.flag().is_mapped() {
            bottomhash.update_dict(
                &(&read.start() + 1),
                read.query_len() as i32,
                &get_umi(&read, separator),
                &read,
            );
        }
    }

    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(&mut self, input_file: BamReader<File>, mut bottomhash: BottomHashMap) {
        let read_puller = match self.group_by_length {
            false => ChunkProcessor::pull_read,
            true => ChunkProcessor::pull_read_w_length,
        };

        let mut counter = 0;
        for r in input_file {
            let read = &r.unwrap();
            read_puller(self, read, &mut bottomhash, self.separator);
            counter += 1;
            if counter % 100_000 == 0 {
                print! {"\rRead in {counter} reads" }
            }
        }
        print! {"\r Grouping {counter} reads...\n"}

        Self::group_reads(self, &mut bottomhash);
    }
}
