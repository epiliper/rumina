use crate::bam_io::bam_reader::WindowedBamReader;
use crate::bottomhash;
use crate::progbars::make_windowbar;
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

// this struct serves to retrieve reads from either indexed or unindexed BAM files, batch them, and
// organize them in the bottomhash data structure for downtream UMI-based operations.
//
pub fn init_processor(
    input_file: String,
    output_file: String,
    grouping_method: GroupingMethod,
    threads: usize,
    strict_threads: bool,
    split_window: Option<i64>,
    group_by_length: bool,
    only_group: bool,
    singletons: bool,
    r1_only: bool,
    min_maxes: Arc<Mutex<GroupReport>>,
    seed: u64,
) -> (
    WindowedBamReader,
    Option<IndexedReader>,
    Writer,
    ChunkProcessor,
) {
    let io_threads = match strict_threads {
        true => threads,
        false => num_cpus::get(),
    };

    // let (header, bam_reader) = make_bam_reader(&input_file, io_threads);
    let bam_reader = WindowedBamReader::new(&input_file, io_threads, split_window);
    let bam_writer = make_bam_writer(&output_file, bam_reader.raw_header.clone(), io_threads);

    let read_handler = ChunkProcessor {
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

    let mate_reader = match r1_only {
        true => Some(make_bam_reader(&input_file, io_threads).1),
        false => None,
    };

    (bam_reader, mate_reader, bam_writer, read_handler)
}

// for every position, group, and process UMIs. output remaining UMIs to write list
pub fn process_chunks(
    chunk_processor: &mut ChunkProcessor,
    mut reader: WindowedBamReader,
    mut other_reader: Option<IndexedReader>,
    separator: &String,
    mut bam_writer: Writer,
) {
    let (mut pos, mut key): (i64, ReadKey);
    let mut outreads: Vec<Record> = Vec::with_capacity(1_000_000);

    let multiprog = MultiProgress::new();

    while reader.next_reference() {
        let mut window_bar = make_windowbar(reader.windows.len() as u64);
        window_bar.set_prefix("WINDOW");
        window_bar = multiprog.add(window_bar);

        while reader.next_window() {
            let mut bottomhash = bottomhash::BottomHashMap {
                read_dict: IndexMap::with_capacity(500),
            };

            let mut window_records = 0;

            for record in reader.window_records() {
                if !record.is_first_in_template() && chunk_processor.r1_only {
                    continue;
                }

                (pos, key) = get_read_pos_key(chunk_processor.group_by_length, &record);
                chunk_processor.pull_read(record, pos, key, &mut bottomhash, &separator);
                window_records += 1;
            }

            info!("{} reads pulled from window", window_records);
            window_bar.set_message(format!("{window_records} reads in window"));
            window_bar.inc(1);

            outreads.extend(ChunkProcessor::group_reads(
                chunk_processor,
                &mut bottomhash,
                &multiprog,
                separator,
            ));

            chunk_processor.write_reads(
                &mut outreads,
                &mut bam_writer,
                &mut other_reader,
                reader.cur_ref,
                &reader.cur_window,
            );

            info!("Processed {} total reads...", chunk_processor.read_counter);
        }
        multiprog.remove(&window_bar);
        window_bar.finish();
    }
}
