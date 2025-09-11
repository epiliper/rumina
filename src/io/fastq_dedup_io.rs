use crate::cli::DedupArgs;
use crate::io::fastqio::{FastqInput, WritesFastqRecords};
use crate::io::FileIO;
use crate::record::FastqRecord;
use anyhow::Error;

use std::path::Path;

use super::fastq_writer::create_writer_from_file;

pub struct FastqIO {
    pub reader: Option<FastqInput>,
    pub _mate_reader: Option<FastqInput>,
    pub writer: Box<dyn WritesFastqRecords>,
    pub _separator: String,
}

impl FastqIO {
    pub fn new(
        infile_name: &str,
        outfile_name: &str,
        retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
        _separator: String,
    ) -> Result<Self, Error> {
        let num_threads = if strict_threads {
            num_threads
        } else {
            num_cpus::get()
        };
        let reader = Some(FastqInput::from_file(Path::new(infile_name))?);

        let _mate_reader = match retrieve_r2s {
            true => Some(FastqInput::from_file(Path::new(infile_name))?),
            false => None,
        };

        let writer = create_writer_from_file(Path::new(outfile_name), num_threads)?;

        Ok(Self {
            reader,
            _mate_reader,
            writer,
            _separator,
        })
    }

    pub fn init_from_args(
        args: &DedupArgs,
        infile_path: &str,
        outfile_path: &str,
    ) -> Result<Self, Error> {
        Self::new(
            infile_path,
            outfile_path,
            args.paired,
            args.threads,
            args.strict_threads,
            args.separator.clone(),
        )
    }
}

impl FileIO<FastqRecord> for FastqIO {
    fn write_reads(&mut self, outreads: &mut Vec<FastqRecord>) {
        for read in outreads.drain(..) {
            self.writer.write_record(read).unwrap()
        }
    }
}
