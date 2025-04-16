use crate::args::Args;
use crate::bam_io::file_io::FileIO;
use crate::read_store::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::{FastqRecord, Record};
use bio::io::fastq::{Reader, Writer};
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};

pub type FastqReader = Reader<BufReader<File>>;
pub type FastqWriter = Writer<BufWriter<File>>;

pub struct FastqIO {
    pub reader: Option<FastqReader>,
    pub _mate_reader: Option<FastqReader>,
    pub writer: FastqWriter,
    pub _num_threads: usize,
    pub separator: String,
}

impl FastqIO {
    pub fn new(
        infile_name: &String,
        outfile_name: &String,
        retrieve_r2s: bool,
        num_threads: usize,
        strict_threads: bool,
        separator: String,
    ) -> Self {
        let _num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };

        let _mate_reader = match retrieve_r2s {
            true => {
                Some(FastqReader::from_file(infile_name).expect("Failed to create mate reader!"))
            }
            false => None,
        };

        let reader =
            Some(Reader::from_file(infile_name).expect("Failed to create reader for file!"));
        let writer = Writer::new(std::io::BufWriter::new(
            File::create(outfile_name).expect("Failed to create output file!"),
        ));

        Self {
            reader,
            _mate_reader,
            writer,
            _num_threads,
            separator,
        }
    }

    pub fn init_from_args(args: &Args, infile_path: &String, outfile_path: &String) -> Self {
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
