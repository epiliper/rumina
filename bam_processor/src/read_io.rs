use crate::bottomhash::{self, BottomHashMap};
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::progbars::*;
use crate::readkey::ReadKey;
use crate::report::{BarcodeTracker, StaticUMI};
use crate::GroupReport;
use crate::GroupingMethod;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Header, IndexedReader, Read, Record, Writer};
use std::collections::HashMap;
use std::sync::Arc;

const WINDOW_CHUNK_SIZE: usize = 3; // number of coord windows processed at once

pub fn get_umi<'b>(record: &'b Record, separator: &String) -> &'b str {
    unsafe {
        std::str::from_utf8_unchecked(record.qname())
            .rsplit_once(separator)
            .expect("ERROR: failed to get UMI from read QNAME. Check --separator. Exiting.")
            .1
    }
}

pub fn get_umi_static<'c>(raw_umi: &'c str) -> StaticUMI {
    let mut umi = StaticUMI::new();
    umi.extend(raw_umi.as_bytes().into_iter().map(|b| *b));
    umi
}

pub fn get_windows(
    window_size: Option<i64>,
    bam: &IndexedReader,
    reference_id: u32,
) -> Vec<[i64; 2]> {
    if let Some(window_size) = window_size {
        let max_pos = bam
            .header()
            .target_len(reference_id)
            .expect("Invalid reference") as i64;

        let mut ranges: Vec<[i64; 2]> = Vec::new();

        let mut j;
        let mut i = 0;

        while i <= max_pos - window_size {
            j = i;
            i = j + window_size;

            ranges.push([j, i]);
        }

        ranges.push([i, max_pos + 1]);

        return ranges;
    } else {
        return vec![[0, i64::MAX]];
    }
}

pub fn make_bam_writer(output_file: &String, header: Header, num_threads: usize) -> Writer {
    let mut bam_writer =
        Writer::from_path(output_file, &header, rust_htslib::bam::Format::Bam).unwrap();
    bam_writer.set_threads(num_threads).unwrap();
    bam_writer
}

pub fn make_bam_reader(input_file: &String, num_threads: usize) -> (Header, IndexedReader) {
    let mut bam_reader = IndexedReader::from_path(input_file).unwrap();
    bam_reader.set_threads(num_threads).unwrap();
    let header = Header::from_template(bam_reader.header());

    (header, bam_reader)
}

// this struct serves to retrieve reads from either indexed or unindexed BAM files, batch them, and
// organize them in the bottomhash data structure for downtream UMI-based operations.
pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub read_counter: i64,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub only_group: bool,
    pub singletons: bool,
    pub track_barcodes: bool,
    pub barcode_tracker: Arc<Mutex<BarcodeTracker>>,
}

impl<'a> ChunkProcessor<'a> {
    pub fn get_read_pos_key(&self, read: &Record) -> (i64, ReadKey) {
        let mut pos;
        let key: ReadKey;

        if read.is_reverse() {
            // set end pos as start to group with forward-reads covering same region
            pos = read.reference_end();
            pos += read.cigar().trailing_softclips(); // pad with right-side soft clip
            key = ReadKey {
                length: read.seq_len() * self.group_by_length as usize,
                reverse: true,
                chr: read.tid() as usize,
            };
            (pos, key)
        } else {
            pos = read.reference_start();
            pos -= read.cigar().leading_softclips(); // pad with left-side soft clip

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
    pub fn group_reads(
        &mut self,
        bottomhash: &mut BottomHashMap,
        multiprog: &MultiProgress,
    ) -> Vec<Record> {
        let grouping_method = Arc::new(&self.grouping_method);

        let mut coord_bar = make_coordbar(bottomhash.read_dict.len() as u64);
        coord_bar = multiprog.add(coord_bar);

        coord_bar.set_prefix("REFERENCE COORDINATE");

        let outreads = Arc::new(Mutex::new(Vec::new()));

        bottomhash.read_dict.par_drain(..).for_each(|position| {
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
                let mut counts: HashMap<&str, i32> = HashMap::with_capacity(umis_reads.len());

                // get number of reads for each raw UMI
                let mut num_umis = 0;

                for umi in &umis {
                    counts.entry(umi.as_str()).or_insert(umis_reads[umi].count);
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

                let tagged_reads = group_handler.tag_records(
                    groupies,
                    umis_reads,
                    Arc::clone(&self.barcode_tracker),
                );

                drop(grouper);

                // update grouping report
                if let Some(tagged_reads) = tagged_reads {
                    let mut out = outreads.lock();
                    out.extend(tagged_reads.1);
                    drop(out);

                    if let Some(group_report) = tagged_reads.0 {
                        let mut min_max = self.min_max.lock();
                        min_max.update(group_report, num_umis);
                        drop(min_max)
                    }
                }
            }
            self.barcode_tracker.lock().write_to_report_file();

            coord_bar.inc(1);
        });

        coord_bar.finish_and_clear();
        Arc::try_unwrap(outreads)
            .expect("Unable to dereference tagged reads!")
            .into_inner()
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

    pub fn write_reads(&mut self, mut outreads: Vec<Record>, bam_writer: &mut Writer) {
        outreads
            .drain(..)
            .for_each(|read| bam_writer.write(&read).unwrap());
    }

    // for every position, group, and process UMIs. output remaining UMIs to write list
    pub fn process_chunks(
        &mut self,
        mut reader: IndexedReader,
        window_size: Option<i64>,
        mut bam_writer: Writer,
    ) {
        let mut pos;
        let mut key;

        let ref_count = reader.header().clone().target_count();

        let mut read_bar = make_readbar();

        for tid in 0..ref_count {
            let windows = get_windows(window_size, &reader, tid);

            let multiprog = MultiProgress::new();
            read_bar = multiprog.add(read_bar);

            let mut window_bar = make_windowbar(windows.len() as u64);
            window_bar = multiprog.add(window_bar);

            window_bar.set_prefix("WINDOW");

            for window_chunk in windows.chunks(WINDOW_CHUNK_SIZE) {
                let mut bottomhash = bottomhash::BottomHashMap {
                    read_dict: IndexMap::with_capacity(500),
                };

                let mut window_reads = 0;
                for window in window_chunk {
                    let start = window[0];
                    let end = window[1];

                    reader
                        .fetch((tid, window[0], window[1]))
                        .expect("Error: invalid window value supplied!");

                    for read in reader.records().map(|read| read.unwrap()) {
                        if read.is_reverse() {
                            // reverse-mapping reads
                            if read.reference_end() <= end && read.reference_end() >= start {
                                (pos, key) = self.get_read_pos_key(&read);
                                self.pull_read(read, pos, key, &mut bottomhash, self.separator);
                                window_reads += 1;
                            }
                            // forward-mapping reads
                        } else if read.reference_start() < end && read.reference_start() >= start {
                            (pos, key) = self.get_read_pos_key(&read);
                            self.pull_read(read, pos, key, &mut bottomhash, self.separator);
                            window_reads += 1;
                        }
                    }
                }
                window_bar.set_message(format!("{window_reads} reads in window"));
                let outreads = Self::group_reads(self, &mut bottomhash, &multiprog);
                bottomhash.read_dict.clear();
                self.write_reads(outreads, &mut bam_writer);
                window_bar.inc(WINDOW_CHUNK_SIZE as u64);
                read_bar.set_message(format!("Processed {} total reads...", self.read_counter));
            }
            window_bar.finish();
        }
        read_bar.finish();
    }
}
