#![allow(unused_variables)]
use anyhow::{Context, Error};
use bio::io::fastq::{Reader, Record, Records, Writer};
use flate2::read::GzDecoder;
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
};

use gzp::{deflate::Gzip, BgzfSyncReader, ZBuilder, ZWriter};

use crate::record::FastqRecord;

/// Used to impose order on reads otherwise unordered. Used right now in channels where reads are inbound from multiple threads.
pub struct IntakeOrdered {
    pub pair: ReadPair,
    pub order: u32,
}

pub struct ReadPair {
    pub r1: Option<Record>,
    pub r2: Option<Record>,
}

/// An interface for iterating over bam records from a gzipped or plaintext fastq.
pub trait FastqRecordIterator {
    fn next_record(&mut self) -> Option<Result<FastqRecord, Error>>;
}

/// for reading plaintext fastq files
impl FastqRecordIterator for Records<BufReader<File>> {
    fn next_record(&mut self) -> Option<Result<FastqRecord, Error>> {
        self.next().map(|e| e.map_err(anyhow::Error::msg))
    }
}

/// For reading gzipped-fastq files
impl FastqRecordIterator for Records<BufReader<GzDecoder<File>>> {
    fn next_record(&mut self) -> Option<Result<FastqRecord, Error>> {
        self.next().map(|e| e.map_err(anyhow::Error::msg))
    }
}

pub type CompressedFastqWriter = Box<dyn ZWriter>;
pub type DecompressedFastqWriter = Writer<File>;

/// an interface for writing compressed and plaintext fastq records to file.
pub trait WritesFastqRecords {
    fn write_record(&mut self, r: FastqRecord) -> Result<(), Error>;
}

impl WritesFastqRecords for CompressedFastqWriter {
    /// A trait implementation that basically tells a gzp writer how to write FastQ records to file.
    /// this is basically stolen from [bio::io::Fastq::Reader::write()].
    fn write_record(&mut self, r: Record) -> Result<(), Error> {
        self.write_all(b"@")?;
        self.write_all(r.id().as_bytes())?;
        if let Some(desc) = r.desc() {
            self.write_all(b" ")?;
            self.write_all(desc.as_bytes())?;
        }
        self.write_all(b"\n")?;
        self.write_all(r.seq())?;
        self.write_all(b"\n+\n")?;
        self.write_all(r.qual())?;
        self.write_all(b"\n")?;

        Ok(())
    }
}

impl WritesFastqRecords for DecompressedFastqWriter {
    fn write_record(&mut self, r: Record) -> Result<(), Error> {
        self.write(r.id(), r.desc(), r.seq(), r.qual())
            .map_err(Error::msg)
    }
}

fn fastq_create_reader_compressed(
    infile: PathBuf,
) -> Result<Records<BufReader<GzDecoder<File>>>, Error> {
    let fhandle = File::open(infile)?;
    let decoder = GzDecoder::new(fhandle);
    let reader = Reader::new(decoder);

    Ok(reader.records())
}

fn fastq_create_reader_decompressed(infile: PathBuf) -> Result<Records<BufReader<File>>, Error> {
    let fhandle = File::open(infile)?;
    let reader = Reader::new(fhandle);

    Ok(reader.records())
}

pub fn fastq_create_writer_compressed(
    outfile: &Path,
    threads: usize,
) -> Result<CompressedFastqWriter, Error> {
    let fhandle = File::create(outfile)?;
    let bufwriter = BufWriter::new(fhandle);

    let writer = ZBuilder::<Gzip, _>::new()
        .num_threads(4)
        .from_writer(bufwriter);
    Ok(writer)
}

pub fn fastq_create_writer_decompressed(
    outfile: &Path,
    _threads: usize,
) -> Result<DecompressedFastqWriter, Error> {
    let fhandle = File::create(outfile)?;
    let bufwriter = BufWriter::new(fhandle);
    let writer = Writer::from_bufwriter(bufwriter);

    Ok(writer)
}

/// A struct representing a record stream from a single fastq file.
///
/// It holds an interator over records of the file, generated from either a plaintext or gzip reader.
pub struct FastqInput {
    _infile: PathBuf,
    pub records: Box<dyn FastqRecordIterator>,
}

impl FastqInput {
    /// Create an input record iterator, detecting whether or not records need to be decompressed.
    pub fn from_file(infile: &Path) -> Result<Self, Error> {
        // Create reader
        let is_gzip = match infile
            .to_str()
            .context("failed to convert file name to string")?
            .rsplit_once(".")
        {
            None => anyhow::bail!("Unable to detect file format!"),
            Some((_pre, end)) => end == "gz",
        };

        match is_gzip {
            true => {
                let records = fastq_create_reader_compressed(infile.to_path_buf())?;
                Ok(Self {
                    records: Box::new(records),
                    _infile: infile.to_path_buf(),
                })
            }
            false => {
                let records = fastq_create_reader_decompressed(infile.to_path_buf())?;
                Ok(Self {
                    records: Box::new(records),
                    _infile: infile.to_path_buf(),
                })
            }
        }
    }
}
