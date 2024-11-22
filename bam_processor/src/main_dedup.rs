use crate::bottomhash;
use crate::progbars::*;
use crate::report::{init_barcode_writer, BarcodeWriter};
use crate::utils::{get_read_pos_key, get_windows, make_bam_reader, make_bam_writer};
use crate::window_processor::*;
use crate::BarcodeTracker;
use crate::GroupReport;
use crate::GroupingMethod;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use parking_lot::Mutex;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{IndexedReader, Read, Writer};
use std::sync::Arc;

const WINDOW_CHUNK_SIZE: usize = 3; // number of coord windows processed at once

// this struct serves to retrieve reads from either indexed or unindexed BAM files, batch them, and
// organize them in the bottomhash data structure for downtream UMI-based operations.
//
pub fn init_processor(
    input_file: String,
    output_file: String,
    grouping_method: GroupingMethod,
    threads: usize,
    split_window: Option<i64>,
    group_by_length: bool,
    only_group: bool,
    singletons: bool,
    track_barcodes: Option<&String>,
    min_maxes: Arc<Mutex<GroupReport>>,
    barcode_tracker: Arc<Mutex<BarcodeTracker>>,
    seed: u64,
) -> (IndexedReader, Writer, ChunkProcessor) {
    let (header, bam_reader) = make_bam_reader(&input_file, threads);
    let bam_writer = make_bam_writer(&output_file, header, threads);

    let read_handler = ChunkProcessor {
        min_max: Arc::clone(&min_maxes),
        grouping_method,
        group_by_length,
        seed,
        split_window,
        only_group,
        singletons,
        read_counter: 0,
        track_barcodes: track_barcodes.cloned(),
        barcode_tracker: Arc::clone(&barcode_tracker),
    };

    (bam_reader, bam_writer, read_handler)
}

// for every position, group, and process UMIs. output remaining UMIs to write list
pub fn process_chunks(
    chunk_processor: &mut ChunkProcessor,
    mut reader: IndexedReader,
    separator: &String,
    mut bam_writer: Writer,
) {
    let mut pos;
    let mut key;

    let (bc_sender, bc_writer): BarcodeWriter;

    if let Some(outfile) = &chunk_processor.track_barcodes {
        (bc_sender, bc_writer) = init_barcode_writer(outfile)
    } else {
        (bc_sender, bc_writer) = (None, None)
    }

    let ref_count = reader.header().clone().target_count();

    let mut read_bar = make_readbar();

    for tid in 0..ref_count {
        let windows = get_windows(chunk_processor.split_window, &reader, tid);

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
                            (pos, key) = get_read_pos_key(chunk_processor.group_by_length, &read);
                            chunk_processor.pull_read(read, pos, key, &mut bottomhash, &separator);
                            window_reads += 1;
                        }
                        // forward-mapping reads
                    } else if read.reference_start() < end && read.reference_start() >= start {
                        (pos, key) = get_read_pos_key(chunk_processor.group_by_length, &read);
                        chunk_processor.pull_read(read, pos, key, &mut bottomhash, &separator);
                        window_reads += 1;
                    }
                }
            }
            window_bar.set_message(format!("{window_reads} reads in window"));
            let outreads = ChunkProcessor::group_reads(
                chunk_processor,
                &mut bottomhash,
                &bc_sender,
                &multiprog,
                separator,
            );
            bottomhash.read_dict.clear();
            chunk_processor.write_reads(outreads, &mut bam_writer);
            window_bar.inc(WINDOW_CHUNK_SIZE as u64);
            read_bar.set_message(format!(
                "Processed {} total reads...",
                chunk_processor.read_counter
            ));
        }
        window_bar.finish();
    }
    if chunk_processor.track_barcodes.is_some() {
        drop(bc_sender.unwrap());
        let _ = bc_writer.unwrap().join();
    }
    read_bar.finish();
}
