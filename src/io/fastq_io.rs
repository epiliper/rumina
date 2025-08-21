use crate::args::Args;
use crate::io::FileIO;
use crate::record::FastqRecord;
use anyhow::Error;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use seq_io::fastq::Reader;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

pub type FastqReaderGz = Reader<BufReader<GzDecoder<File>>>;
pub type FastqWriter = GzEncoder<BufWriter<File>>;

pub struct FastqIO {
    pub reader: Option<FastqReaderGz>,
    pub _mate_reader: Option<FastqReaderGz>,
    pub writer: FastqWriter,
    pub _num_threads: usize,
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
        let _num_threads = match strict_threads {
            true => num_threads,
            false => num_cpus::get(),
        };

        let inf = File::open(infile_name)?;
        let decoder = GzDecoder::new(inf);
        let buf_reader = BufReader::with_capacity(8192, decoder);

        let _mate_reader = match retrieve_r2s {
            true => {
                let inf = File::open(infile_name)?;
                let decoder = GzDecoder::new(inf);
                let buf_reader = BufReader::with_capacity(8192, decoder);
                Some(Reader::new(buf_reader))
            }
            false => None,
        };

        let reader = Some(Reader::new(buf_reader));

        let outf = File::create(outfile_name)?;
        let buf_writer = BufWriter::with_capacity(8192, outf);
        let writer = GzEncoder::new(buf_writer, Compression::default());

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

    pub fn write_record(&mut self, record: FastqRecord) -> Result<(), Error> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(&record.head)?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(&record.seq)?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(&record.qual)?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }
}

impl FileIO<FastqRecord> for FastqIO {
    fn write_reads(&mut self, outreads: &mut Vec<FastqRecord>) {
        for read in outreads.drain(..) {
            self.write_record(read).unwrap()
        }
    }
}
