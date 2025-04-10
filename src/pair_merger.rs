use crate::merge::handle_dupes;
use crate::merge_report::MergeReport;
use crate::read_store::pair_bundles::*;
use crate::realign::init_remapper;
use crate::utils::{get_windows, make_bam_reader, make_bam_writer};
use crossbeam::channel::{unbounded, Receiver, Sender};
use indexmap::IndexMap;
use log::info;
use rust_htslib::bam::{Read, Writer};
use std::{thread, thread::JoinHandle};

use rust_htslib::bam::Record;

pub fn spawn_writer_thread(mut bam_writer: Writer, r: Receiver<Record>) -> JoinHandle<i32> {
    thread::spawn(move || {
        let mut buffer: Vec<Record> = Vec::with_capacity(1_000_000);
        let mut num_writes = 0;

        loop {
            match r.recv() {
                Ok(read) => buffer.push(read),
                Err(_) => {
                    buffer.sort_by(|ra, rb| ra.pos().cmp(&rb.pos()));
                    for read in buffer.drain(..) {
                        bam_writer.write(&read).expect("unable to write read");
                        num_writes += 1;
                    }
                    buffer.clear();

                    return num_writes;
                }
            }
        }
    })
}

#[derive(Debug)]
pub struct PairMerger {
    pub ref_fasta: String,
    pub min_overlap_bp: i64,
    pub threads: usize,
    pub infile: String,
    pub outfile: String,
    pub split_window: Option<i64>,
}

impl PairMerger {
    pub fn merge_windows(&mut self) -> MergeReport {
        let mut merge_report = MergeReport::new();

        let (header, mut reader) = make_bam_reader(&self.infile, self.threads);
        let (mapper, ref_fasta) = init_remapper(&self.ref_fasta);
        let mut num_writes: i32 = 0;

        let ref_count = reader.header().clone().target_count();
        let mut read_count = 0;
        for tid in 0..ref_count {
            let max_pos = reader.header().target_len(tid).unwrap() as i64;
            let windows = get_windows(self.split_window, max_pos);
            reader.fetch((tid, 0, u32::MAX)).unwrap();
            let mut next_window_reads: Vec<Record> = Vec::with_capacity(100);

            for window_chunk in windows.chunks(3) {
                let mut bundles = PairBundles {
                    read_dict: IndexMap::new(),
                };

                let writer = make_bam_writer(&self.outfile, header.clone(), self.threads);
                let (s, r): (Sender<Record>, Receiver<Record>) = unbounded();
                let writer_handle = spawn_writer_thread(writer, r);

                for window in window_chunk {
                    let start = window.start;
                    let end = window.end;

                    info!("Ref: {}, Start: {}, End: {}", tid, start, end);

                    next_window_reads.drain(..).for_each(|record| {
                        bundles.update_dict(record);
                        read_count += 1;
                    });

                    for read in reader.records().flatten() {
                        if read.pos() >= end {
                            next_window_reads.push(read);
                            break;
                        } else if read.pos() < end && read.pos() >= start {
                            bundles.update_dict(read);
                            read_count += 1;
                        }
                    }
                }

                let merge_results = handle_dupes(
                    &mut bundles.read_dict,
                    mapper.clone(),
                    &ref_fasta,
                    self.min_overlap_bp as usize,
                    s.clone(),
                );

                for res in merge_results {
                    merge_report.count(res);
                }

                drop(s);
                num_writes += writer_handle.join().expect("Writer thread panicked");
            }
        }

        merge_report.num_inreads = read_count;
        merge_report.num_outreads = num_writes;
        info!("{:?}", merge_report);

        merge_report
    }
}
