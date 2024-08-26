use crate::bottomhash::BottomHashMap;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::GroupingMethod;
use indicatif::ProgressBar;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, Record};
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::HashMap;
use std::sync::Arc;

fn get_umi(record: &Record, separator: &String) -> String {
    unsafe {
        std::str::from_utf8_unchecked(record.qname())
            .rsplit_once(separator)
            .unwrap()
            .1
            .to_string()
    }
}

pub fn get_read_pos(read: &Record) -> Option<i64> {
    let mut pos;

    if !read.cigar().is_empty() {
        pos = read.reference_end();

        if read.is_reverse() {
            // set end pos as start to group with forward-reads covering same region
            pos += read.cigar().leading_softclips(); // pad with right-side soft clip
            return Some(pos);
        } else {
            pos = read.pos();

            pos -= read.cigar().trailing_softclips(); // pad with left-side soft clip

            return Some(pos);
        };
    }

    return None;
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
    pub num_reads_input_file: i64,
    pub num_reads_output_file: i64,
}

pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub read_counter: i64,
    pub reads_to_output: Arc<Mutex<Vec<Record>>>,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub only_group: bool,
    pub singletons: bool,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam

    pub fn group_reads(&mut self, bottomhash: &mut BottomHashMap, drain_end: usize) {
        let grouping_method = Arc::new(&self.grouping_method);

        let current_len = bottomhash.bottom_dict.len();

        let range = match drain_end {
            0 => 0..current_len,
            _ => 0..std::cmp::min(drain_end, current_len),
        };

        bottomhash
            .bottom_dict
            .par_drain(range)
            .for_each(|position| {
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

                    let grouper = Grouper { umis: &umis };
                    let mut counts: HashMap<&String, i32> =
                        HashMap::with_capacity(umis_reads.len());

                    // get number of reads for each raw UMI
                    let mut num_umis = 0;
                    for umi in &umis {
                        counts.entry(umi).or_insert(umis_reads[umi].count);
                        num_umis += 1;
                    }

                    let mut group_handler = GroupHandler {
                        seed: self.seed + position.0 as u64, // make seed unique per position
                        group_only: self.only_group,
                        singletons: self.singletons,
                    };

                    // perform UMI clustering per the method specified
                    let groupies = grouper.cluster(counts, Arc::clone(&grouping_method));

                    let tagged_reads = group_handler.tag_records(groupies, umis_reads);
                    let mut min_max = self.min_max.lock();

                    // update the groups with mininum and maximum observed reads
                    match tagged_reads {
                        Some(tagged_reads) => {
                            self.reads_to_output.lock().extend(tagged_reads.1);
                            // .expect("Read channel not recieving!");

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

                                    // record the number of reads to be written
                                    min_max.num_reads_output_file += x.num_reads_output_file;
                                }
                                _ => (),
                            }
                        }
                        None => (),
                    }
                }
            });
    }

    // organize reads in bottomhash based on position
    pub fn pull_read(
        &mut self,
        read: &Record,
        pos: i32,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        bottomhash.update_dict(pos, 0, &get_umi(&read, separator), &read);
        self.read_counter += 1;
    }

    pub fn pull_read_w_length(
        &mut self,
        read: &Record,
        pos: i32,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        bottomhash.update_dict(read.seq_len() as i32, pos, &get_umi(&read, separator), read);

        self.read_counter += 1;
    }

    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(&mut self, mut input_file: IndexedReader, mut bottomhash: BottomHashMap) {
        let progressbar = ProgressBar::new(bottomhash.bottom_dict.keys().len().try_into().unwrap());

        let read_puller = match self.group_by_length {
            false => ChunkProcessor::pull_read,
            true => ChunkProcessor::pull_read_w_length,
        };

        let mut pos;

        let mut read: Record;

        // todo: make this neater
        input_file.fetch(FetchDefinition::All).unwrap();
        for r in input_file.records() {
            read = r.unwrap();

            if read.is_unmapped() {
                continue;
            }

            // get adjusted pos for bottomhash
            // use start to keep track of position in BAM file.
            pos = get_read_pos(&read).expect("ERROR: invalid CIGAR string!!!");

            read_puller(self, &read, pos as i32, &mut bottomhash, self.separator);

            if self.read_counter % 100_000 == 0 {
                print! {"\rRead in {} reads", self.read_counter}
            }
        }
        print! {"\r Grouping {} reads...\n", self.read_counter}

        println!("processing remaining reads...");
        Self::group_reads(self, &mut bottomhash, 0);
        progressbar.finish();
    }
}
