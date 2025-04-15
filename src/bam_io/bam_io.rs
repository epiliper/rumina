use crate::args::Args;
use crate::bam_io::bam_reader::WindowedBamReader;
use crate::utils::{make_bam_reader, make_bam_writer};
use rust_htslib::bam::{IndexedReader, Writer};

pub struct BamIO {
    pub windowed_reader: crate::bam_io::bam_reader::WindowedBamReader,
    pub mate_reader: Option<IndexedReader>,
    pub writer: Writer,
    pub num_threads: usize,
    pub window_size: Option<i64>,
}

impl BamIO {
    pub fn new(
        infile_name: &String,
        outfile_name: &String,
        retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
        window_size: Option<i64>,
    ) -> Self {
        let num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };
        let windowed_reader = WindowedBamReader::new(infile_name, num_threads, window_size);
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
            window_size,
        }
    }

    pub fn init_from_args(args: &Args, infile_path: &String, outfile_path: &String) -> Self {
        Self::new(
            infile_path,
            outfile_path,
            args.r1_only,
            args.threads,
            args.strict_threads,
            args.split_window,
        )
    }
}
