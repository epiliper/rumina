use crate::bottomhash::ReadsAndCount;
use crate::main_dedup::WINDOW_CHUNK_SIZE;
use crate::merge::handle_dupes;
use crate::merge_report::MergeReport;
use crate::realign::init_remapper;
use crate::utils::{get_windows, make_bam_reader, make_bam_writer};
use crossbeam::channel::{bounded, Receiver, Sender};
use indexmap::IndexMap;
use log::{info, warn};
use rust_htslib::bam::{record::Aux, Read};
use std::thread;

use rust_htslib::bam::Record;
pub struct PairBundles {
    read_dict: IndexMap<String, ReadsAndCount>,
}

const UMI_TAG: &[u8; 2] = b"BX";

impl PairBundles {
    pub fn update_dict(&mut self, umi: String, read: Record) {
        self.read_dict
            .entry(umi)
            .or_insert_with(|| ReadsAndCount {
                reads: Vec::new(),
                count: 0,
            })
            .up(read)
    }
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
        let mut writer = make_bam_writer(&self.outfile, header, self.threads);
        let (mapper, ref_fasta) = init_remapper(&self.ref_fasta);

        let (s, r): (Sender<Record>, Receiver<Record>) = bounded(10000);

        let writer_handle = thread::spawn(move || {
            let mut buffer = Vec::with_capacity(1000000);
            let mut counter: u32 = 0;
            let mut num_writes = 0;

            loop {
                match r.recv() {
                    Ok(read) => {
                        buffer.push(read);
                        counter += 1;
                        if counter == 1000 {
                            for read in &buffer {
                                writer.write(read).expect("Error writing read");
                                num_writes += 1;
                            }
                            buffer.clear();
                            counter = 0;
                        }
                    }
                    Err(_) => {
                        if !buffer.is_empty() {
                            for read in &buffer {
                                writer.write(read).expect("Error writing read");
                                num_writes += 1;
                            }
                        }
                        return num_writes;
                    }
                }
            }
        });

        let ref_count = reader.header().clone().target_count();
        let mut read_count = 0;
        for tid in 0..ref_count {
            let windows = get_windows(self.split_window, &reader, tid);

            for window_chunk in windows.chunks(WINDOW_CHUNK_SIZE) {
                let mut bundles = PairBundles {
                    read_dict: IndexMap::new(),
                };

                for window in window_chunk {
                    reader
                        .fetch((tid, window[0], window[1]))
                        .expect("Error: invalid window value supplied!");

                    for read in reader.records().map(|read| read.unwrap()) {
                        if read.pos() < window[1] && read.pos() >= window[0] {
                            let umi = if let Ok(Aux::String(bx_i)) = read.aux(UMI_TAG) {
                                bx_i
                            } else {
                                warn!("Cannot find UMI for read: {:?}", read);
                                "NULL"
                            };

                            bundles.update_dict(umi.to_string(), read);
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
            }
        }

        drop(s);
        let num_writes = writer_handle.join().expect("Writer thread panicked!");
        merge_report.num_inreads = read_count;
        merge_report.num_outreads = num_writes;
        info!("{:?}", merge_report);

        merge_report
    }
}
