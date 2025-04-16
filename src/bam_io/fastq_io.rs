use crate::args::Args;
use crate::record::FastqRecord;
use bio::io::fastq::{Reader, Writer};
use std::fs::File;
use std::io::{BufReader, BufWriter};

pub type FastqReader = Reader<BufReader<File>>;
pub type FastqWriter = Writer<BufWriter<File>>;

pub struct FastqIO {
    pub reader: FastqReader,
    pub writer: FastqWriter,
    pub num_threads: usize,
}

impl FastqIO {
    pub fn new(
        infile_name: &String,
        outfile_name: &String,
        _retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
    ) -> Self {
        let num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };

        let reader = Reader::from_file(infile_name).expect("Failed to create reader for file!");
        let writer = Writer::new(std::io::BufWriter::new(
            File::create(outfile_name).expect("Failed to create output file!"),
        ));

        Self {
            reader,
            writer,
            num_threads,
        }
    }

    pub fn init_from_args(args: &Args, infile_path: &String, outfile_path: &String) -> Self {
        Self::new(
            infile_path,
            outfile_path,
            args.r1_only,
            args.threads,
            args.strict_threads,
        )
    }
}
