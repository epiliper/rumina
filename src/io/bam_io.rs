use crate::args::Args;
use crate::io::{FileIO, WindowedBamReader};
use crate::record::BamRecord;
use crate::utils::{make_bam_reader, make_bam_writer};
use indexmap::IndexSet;
use log::info;
use rayon::prelude::ParallelSliceMut;
use rust_htslib::bam::{IndexedReader, Read, Writer};
use std::cmp;

pub struct BamIO {
    pub windowed_reader: WindowedBamReader,
    pub mate_reader: Option<IndexedReader>,
    pub writer: Writer,
    pub num_threads: usize,
    pub _window_size: Option<i64>,
    pub _separator: String,
}

impl BamIO {
    pub fn new(
        infile_name: &String,
        outfile_name: &String,
        retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
        _window_size: Option<i64>,
        _separator: String,
    ) -> Self {
        let num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };
        let windowed_reader = WindowedBamReader::new(infile_name, num_threads, _window_size);
        let writer = make_bam_writer(
            outfile_name,
            windowed_reader.raw_header.clone(),
            num_threads,
        );
        let mate_reader = match retrieve_r2s {
            true => Some(make_bam_reader(infile_name, num_threads).1),
            false => None,
        };

        Self {
            windowed_reader,
            mate_reader,
            writer,
            num_threads,
            _window_size,
            _separator,
        }
    }

    pub fn init_from_args(args: &Args, infile_path: &String, outfile_path: &String) -> Self {
        Self::new(
            infile_path,
            outfile_path,
            args.paired,
            args.threads,
            args.strict_threads,
            args.split_window,
            args.separator.clone(),
        )
    }

    pub fn retrieve_r2s(&mut self, ids: IndexSet<&[u8]>) -> Option<Vec<BamRecord>> {
        let ref_len = self
            .windowed_reader
            .meta_header
            .target_len(self.windowed_reader.cur_ref)
            .unwrap_or(u64::MAX);

        if let Some(mate_reader) = &mut self.mate_reader {
            mate_reader
                .fetch((
                    self.windowed_reader.cur_ref,
                    cmp::min(0, self.windowed_reader.cur_window.start - 100),
                    cmp::max(self.windowed_reader.cur_window.end + 100, ref_len as i64),
                ))
                .unwrap();

            let mut mates: Vec<BamRecord> = Vec::with_capacity(ids.len());
            for read in mate_reader.records().flatten() {
                if ids.contains(read.qname()) && read.is_last_in_template() {
                    mates.push(read);
                }
            }

            Some(mates)
        } else {
            None
        }
    }
}

impl FileIO<BamRecord> for BamIO {
    fn write_reads(&mut self, outreads: &mut Vec<BamRecord>) {
        let mut count = 0;
        let mut mates: Option<Vec<BamRecord>> = None;

        if !outreads.is_empty() {
            if let Some(ref mut _bam_reader) = self.mate_reader {
                let mut ids_to_pair = IndexSet::with_capacity(outreads.len());

                outreads.iter().for_each(|read| {
                    ids_to_pair.insert(read.qname());
                });

                mates = self.retrieve_r2s(ids_to_pair);
            }

            if let Some(mates) = mates {
                outreads.extend(mates);
            }

            outreads.par_sort_by(|ra, rb| ra.pos().cmp(&rb.pos()));
            outreads.drain(..).for_each(|read| {
                self.writer.write(&read).unwrap();
                count += 1;
            });
        }
        info!("Written {count} reads!")
    }
}
