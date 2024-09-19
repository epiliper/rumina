use crate::bottomhash::BottomHashMap;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::readkey::ReadKey;
use crate::report::StaticUMI;
use crate::GroupReport;
use crate::GroupingMethod;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Header, IndexedReader, Read, Reader, Record, Writer};
use rust_htslib::htslib;
use std::collections::HashMap;
use std::sync::Arc;

pub enum BamReader {
    Indexed(IndexedReader),
    Standard(Reader),
}

pub fn get_umi<'a>(record: &'a Record, separator: &String) -> &'a str {
    unsafe {
        std::str::from_utf8_unchecked(record.qname())
            .rsplit_once(separator)
            .expect("ERROR: failed to get UMI from read QNAME. Check --separator. Exiting.")
            .1
        // .to_string()
    }
}

pub fn get_umi_static<'a>(raw_umi: &'a str) -> StaticUMI {
    let mut umi = StaticUMI::new();
    umi.extend(
        raw_umi.as_bytes().into_iter().map(|b| *b), // .into_iter()
                                                    // .map(|b| *b),
    );

    umi
}

pub fn make_bam_writer(output_file: &String, header: Header, num_threads: usize) -> Writer {
    let mut bam_writer =
        Writer::from_path(output_file, &header, rust_htslib::bam::Format::Bam).unwrap();
    bam_writer.set_threads(num_threads).unwrap();
    bam_writer
}

pub fn make_bam_reader(
    input_file: &String,
    indexed: bool,
    num_threads: usize,
) -> (Header, BamReader) {
    if indexed {
        let mut bam_reader = IndexedReader::from_path(input_file).unwrap();
        bam_reader.set_threads(num_threads).unwrap();
        let header = Header::from_template(bam_reader.header());

        (header, BamReader::Indexed(bam_reader))
    } else {
        let mut bam_reader = Reader::from_path(input_file).unwrap();
        bam_reader.set_threads(num_threads).unwrap();
        let header = Header::from_template(bam_reader.header());

        (header, BamReader::Standard(bam_reader))
    }
}

// this struct serves to retrieve reads from either indexed or unindexed BAM files, batch them, and
// organize them in the bottomhash data structure for downtream UMI-based operations.
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
    pub track_barcodes: bool,
}

impl<'a> ChunkProcessor<'a> {
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
                chr: read.tid() as usize,
            };
            (pos, key)
        } else {
            pos = read.reference_start();
            pos -= read.cigar().trailing_softclips(); // pad with left-side soft clip

            key = ReadKey {
                length: read.seq_len() * self.group_by_length as usize,
                reverse: false,
                chr: read.tid() as usize,
            };
            (pos, key)
        }
    }

    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(&mut self, bottomhash: &mut BottomHashMap, drain_end: usize) {
        let grouping_method = Arc::new(&self.grouping_method);

        let current_len = bottomhash.read_dict.len();

        let range = match drain_end {
            0 => 0..current_len,
            _ => 0..std::cmp::min(drain_end, current_len),
        };

        bottomhash.read_dict.par_drain(range).for_each(|position| {
            for umi in position.1 {
                let mut umis_reads = umi.1;

                // sort UMIs by read count
                // note that this is an unstable sort, so we need to identify read-tied groups
                umis_reads
                    .par_sort_by(|_umi1, count1, _umi2, count2| count2.count.cmp(&count1.count));

                let umis = umis_reads
                    .keys()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();

                let grouper = Grouper { umis: &umis };
                let mut counts: HashMap<&String, i32> = HashMap::with_capacity(umis_reads.len());

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
                    separator: self.separator,
                    track_barcodes: self.track_barcodes,
                };

                // perform UMI clustering per the method specified
                let groupies = grouper.cluster(counts, Arc::clone(&grouping_method));

                let tagged_reads = group_handler.tag_records(groupies, umis_reads);
                let mut min_max = self.min_max.lock();

                // update grouping report
                if let Some(tagged_reads) = tagged_reads {
                    self.reads_to_output.lock().extend(tagged_reads.1);

                    if let Some(group_report) = tagged_reads.0 {
                        min_max.update(group_report, num_umis);
                    }
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
        bottomhash.update_dict(
            pos,
            key.get_key(),
            get_umi(&read, separator).to_string(),
            read,
        );
        self.read_counter += 1;
    }

    pub fn write_reads(&mut self, bam_writer: &mut Writer) {
        self.reads_to_output
            .lock()
            .drain(0..)
            .for_each(|read| bam_writer.write(&read).unwrap())
    }

    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(
        &mut self,
        input_file: BamReader,
        mut bam_writer: Writer,
        mut bottomhash: BottomHashMap,
    ) {
        let mut pos;
        let mut key;

        match input_file {
            BamReader::Standard(mut reader) => {
                for read in reader
                    .records()
                    .map(|x| x.unwrap())
                    .filter(|read| read.flags() & htslib::BAM_FUNMAP as u16 == 0)
                {
                    (pos, key) = self.get_read_pos_key(&read);
                    self.pull_read(read, pos, key, &mut bottomhash, self.separator);

                    if self.read_counter % 100_000 == 0 {
                        print! {"Read in {} reads\r", self.read_counter}
                    }
                }

                println!("Grouping {} reads...", self.read_counter);
                self.group_reads(&mut bottomhash, 0);
                self.write_reads(&mut bam_writer)
            }

            BamReader::Indexed(mut reader) => {
                let ref_count = reader.header().clone().target_count();

                for tid in 0..ref_count {
                    reader.fetch((tid, 0, u32::MAX)).unwrap();

                    for read in reader.records().map(|read| read.unwrap()) {
                        (pos, key) = self.get_read_pos_key(&read);
                        self.pull_read(read, pos, key, &mut bottomhash, self.separator);

                        if self.read_counter % 100_000 == 0 {
                            print! {"Read in {} reads\r", self.read_counter}
                        }
                    }
                }
                println!("Grouping {} reads...", self.read_counter);
                Self::group_reads(self, &mut bottomhash, 0);
                self.write_reads(&mut bam_writer)
            }
        }
    }
}
