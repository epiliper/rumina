use crate::bam_io::{
    bam_reader::WindowedBamReader,
    fastq_io::{FastqReader, FastqWriter},
};
use crate::progbars::make_windowbar;
use crate::read_store::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::{BamRecord, FastqRecord};
use crate::utils::get_read_pos_key;
use crate::window_processor::*;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use log::info;
use rust_htslib::bam::{IndexedReader, Writer};

pub fn process_whole(
    chunk_processor: &mut Processor,
    mut reader: FastqReader,
    separator: &String,
    mut writer: FastqWriter,
) {
    unimplemented!();
}

// for every position, group, and process UMIs. output remaining UMIs to write list
pub fn process_windows(
    chunk_processor: &mut Processor,
    mut reader: WindowedBamReader,
    mut other_reader: Option<IndexedReader>,
    separator: &String,
    mut bam_writer: Writer,
) {
    let (mut pos, mut key): (i64, ReadKey);
    let mut outreads: Vec<BamRecord> = Vec::with_capacity(1_000_000);

    let multiprog = MultiProgress::new();

    while reader.next_reference() {
        let mut window_bar = make_windowbar(reader.windows.len() as u64);
        window_bar.set_prefix("WINDOW");
        window_bar = multiprog.add(window_bar);

        while reader.next_window() {
            let mut bottomhash = BottomHashMap {
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

            outreads.extend(Processor::group_reads(
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
