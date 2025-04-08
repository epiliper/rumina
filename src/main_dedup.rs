use crate::bam_io::bam_reader::WindowedBamReader;
use crate::bottomhash;
use crate::cli::progbars::make_windowbar;
use crate::readkey::ReadKey;
use crate::utils::{get_read_pos_key, make_bam_reader, make_bam_writer};
use crate::window_processor::*;
use crate::GroupReport;
use crate::GroupingMethod;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use log::info;
use parking_lot::Mutex;
use rust_htslib::bam::{IndexedReader, Record, Writer};
use std::sync::Arc;

pub const WINDOW_CHUNK_SIZE: usize = 3; // number of coord windows processed at once
                                        //
                                        // this struct serves to retrieve reads from either indexed or unindexed BAM files, batch them, and
                                        // organize them in the bottomhash data structure for downtream UMI-based operations.
                                        //
                                        //
                                        //
pub trait BamProcessor {
    fn process_bam(&mut self);
}

pub struct BamGroupProcessor {
    pub bam_io: BamIO,
    pub clusterer: Clusterer,
    pub separator: String,
}

impl BamGroupProcessor {
    pub fn new(
        input_file: String,
        output_file: String,
        grouping_method: GroupingMethod,
        threads: usize,
        strict_threads: bool,
        split_window: Option<i64>,
        group_by_length: bool,
        only_group: bool,
        singletons: bool,
        separator: String,
        r1_only: bool,
        min_maxes: Arc<Mutex<GroupReport>>,
        seed: u64,
    ) -> Self {
        let io_threads = match strict_threads {
            true => threads,
            false => num_cpus::get(),
        };

        let bam_io = BamIO::new(&input_file, &output_file, r1_only, io_threads, split_window);

        let clusterer = Clusterer {
            min_max: Arc::clone(&min_maxes),
            grouping_method,
            group_by_length,
            seed,
            split_window,
            only_group,
            singletons,
            read_counter: 0,
            r1_only,
        };

        Self {
            bam_io,
            clusterer,
            separator,
        }
    }
}

impl BamProcessor for BamGroupProcessor {
    fn process_bam(&mut self) {
        let (mut pos, mut key): (i64, ReadKey);
        let mut outreads: Vec<Record> = Vec::with_capacity(1_000_000);

        let multiprog = MultiProgress::new();

        while self.bam_io.windowed_reader.next_reference() {
            let mut window_bar = make_windowbar(self.bam_io.windowed_reader.windows.len() as u64);
            window_bar.set_prefix("WINDOW");
            window_bar = multiprog.add(window_bar);

            while self.bam_io.windowed_reader.next_window() {
                let mut bottomhash = bottomhash::BottomHashMap {
                    read_dict: IndexMap::with_capacity(500),
                };

                let mut window_records = 0;

                for record in self.bam_io.windowed_reader.window_records() {
                    if !record.is_first_in_template() && self.clusterer.r1_only {
                        continue;
                    }

                    (pos, key) = get_read_pos_key(self.clusterer.group_by_length, &record);
                    self.clusterer
                        .pull_read(record, pos, key, &mut bottomhash, &self.separator);
                    window_records += 1;
                }

                info!("{} reads pulled from window", window_records);
                window_bar.set_message(format!("{window_records} reads in window"));
                window_bar.inc(1);

                outreads.extend(self.clusterer.group_reads(
                    &mut bottomhash,
                    &multiprog,
                    &self.separator,
                ));

                self.clusterer.write_reads(
                    &mut outreads,
                    &mut self.bam_io.writer,
                    &mut self.bam_io.mate_reader,
                    self.bam_io.windowed_reader.cur_ref,
                    &self.bam_io.windowed_reader.cur_window,
                );

                info!("Processed {} total reads...", self.clusterer.read_counter);
            }
            multiprog.remove(&window_bar);
            window_bar.finish();
        }
    }
}

pub struct BamIO {
    pub windowed_reader: WindowedBamReader,
    pub mate_reader: Option<IndexedReader>,
    pub writer: Writer,
}

impl BamIO {
    pub fn new(
        input_file: &String,
        output_file: &String,
        r1_only: bool,
        num_threads: usize,
        split_window: Option<i64>,
    ) -> Self {
        let windowed_reader = WindowedBamReader::new(&input_file, num_threads, split_window);
        let writer = make_bam_writer(
            &output_file,
            windowed_reader.raw_header.clone(),
            num_threads,
        );
        let mate_reader = match r1_only {
            true => Some(make_bam_reader(&input_file, num_threads).1),
            false => None,
        };

        Self {
            windowed_reader,
            mate_reader,
            writer,
        }
    }
}
