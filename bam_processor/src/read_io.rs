use crate::bottomhash::BottomHashMap;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::readkey::ReadKey;
use crate::GroupingMethod;
use indicatif::ProgressBar;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, Record};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::htslib;
use std::collections::HashMap;
use std::sync::Arc;

fn get_umi(record: &Record, separator: &String) -> String {
    unsafe {
        std::str::from_utf8_unchecked(record.qname())
            .rsplit_once(separator)
            .expect("ERROR: failed to get UMI from read QNAME. Check --separator. Exiting.")
            .1
            .to_string()
    }
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
    
    pub fn get_read_pos_key(&self, read: &Record) -> (i64, ReadKey) {
        let mut pos;
        let key: ReadKey;

        if read.is_reverse() {
            // set end pos as start to group with forward-reads covering same region
            pos = read.reference_end();
            pos += read.cigar().leading_softclips(); // pad with right-side soft clip
                                                     //
            key = ReadKey {
                length: read.seq_len() * self.group_by_length as usize,
                reverse: true,
            };

            (pos, key)

        } else {
            pos = read.pos();
            pos -= read.cigar().trailing_softclips(); // pad with left-side soft clip

            key = ReadKey {
                length: read.seq_len() * self.group_by_length as usize,
                reverse: false,
            };

            (pos, key)
        }
    }

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
                    umis_reads.par_sort_by(|_umi1, count1, _umi2, count2| {
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
        read: Record,
        pos: i64,
        key: ReadKey,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        bottomhash.update_dict(pos, key.get_key(), get_umi(&read, separator), read);
        self.read_counter += 1;
    }


    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(&mut self, mut input_file: IndexedReader, mut bottomhash: BottomHashMap) {
        let progressbar = ProgressBar::new(bottomhash.bottom_dict.keys().len().try_into().unwrap());

        let mut pos;
        let mut key;
        let mut read;

        // todo: make this neater
        input_file.fetch(FetchDefinition::All).unwrap();

        for r in input_file
            .records()
            .map(|x| x.unwrap())
            .filter(|read| read.flags() & htslib::BAM_FUNMAP as u16 == 0)
        {
            read = r;

            // get adjusted pos for bottomhash
            // use start to keep track of position in BAM file.
            (pos, key) = self.get_read_pos_key(&read);
            self.pull_read(read, pos, key, &mut bottomhash, self.separator);

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
