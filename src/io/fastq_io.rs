use crate::args::Args;
use crate::io::FileIO;
use crate::record::FastqRecord;
use anyhow::{Context, Error};
use bio::io::fastq::{Reader, Writer};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufReader, BufWriter};

pub type FastqReaderGz = Reader<BufReader<GzDecoder<File>>>;
pub type FastqWriter = Writer<BufWriter<File>>;

pub struct FastqIO {
    pub reader: Option<FastqReaderGz>,
    pub _mate_reader: Option<FastqReaderGz>,
    pub writer: FastqWriter,
    pub _num_threads: usize,
    pub _separator: String,
}

impl FastqIO {
    pub fn new(
        infile_name: &String,
        outfile_name: &String,
        retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
        _separator: String,
    ) -> Result<Self, Error> {
        let _num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };

        let _mate_reader = match retrieve_r2s {
            true => Some(Reader::new(GzDecoder::new(
                File::open(infile_name).context("Failed to open mate file")?,
            ))),
            false => None,
        };

        let reader = Some(Reader::new(GzDecoder::new(
            File::open(infile_name).context("Failed to open mate file")?,
        )));

        let writer = Writer::new(std::io::BufWriter::new(
            File::create(outfile_name).context("Failed to create output file!")?,
        ));

        Ok(Self {
            reader,
            _mate_reader,
            writer,
            _num_threads,
            _separator,
        })
    }

    pub fn init_from_args(
        args: &Args,
        infile_path: &String,
        outfile_path: &String,
    ) -> Result<Self, Error> {
        Self::new(
            infile_path,
            outfile_path,
            args.r1_only,
            args.threads,
            args.strict_threads,
            args.separator.clone(),
        )
    }
}

impl FileIO<FastqRecord> for FastqIO {
    fn write_reads(&mut self, outreads: &mut Vec<FastqRecord>) {
        for read in outreads {
            self.writer.write_record(read).unwrap();
        }
    }
}
