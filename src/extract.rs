use crate::{
    io::{
        fastq_extract_io::FastqIO,
        fastqio::{IntakeOrdered, ReadPair},
    },
    ExtractArgs, FastqQualEncoding, PHRED33, PHRED64, SOLEXA,
};

use anyhow::{Context, Error};
use bio::io::fastq::Record;
use crossbeam::channel::{unbounded, Receiver, Sender};

#[derive(Clone, Debug, Eq, PartialEq)]
enum ExtractBaseType {
    Umi,
    CellBarcode,
    Misc,
}

const CELL_BARCODE_CHAR: char = 'C';
const UMI_CHAR: char = 'N';
const MISC_CHAR: char = 'X';

fn parse_detect_string(dstring: &str) -> Result<Vec<ExtractBaseType>, Error> {
    if dstring.chars().into_iter().any(|c| c.is_ascii_digit()) {
        parse_detect_string_alphanumeric(dstring)
    } else {
        parse_detect_string_alphabetic(dstring)
    }
}

fn parse_detect_string_alphabetic(dstring: &str) -> Result<Vec<ExtractBaseType>, Error> {
    let mut out = vec![];
    for c in dstring.chars() {
        match c.to_ascii_uppercase() {
            CELL_BARCODE_CHAR => out.push(ExtractBaseType::CellBarcode),
            UMI_CHAR => out.push(ExtractBaseType::Umi),
            MISC_CHAR => out.push(ExtractBaseType::Misc),
            e => anyhow::bail!(
                "Illegal character specified for pattern {dstring}: {e}; Choose from C, N, or X."
            ),
        }
    }

    Ok(out)
}

fn parse_detect_string_alphanumeric(dstring: &str) -> Result<Vec<ExtractBaseType>, Error> {
    let mut out = vec![];
    let mut intbuf: String = String::with_capacity(4);
    for c in dstring.chars() {
        if c.is_ascii_digit() {
            intbuf.push(c);
        } else {
            if intbuf.is_empty() {
                anyhow::bail!("Pattern character specified before count: {dstring}!\nSpecify count before character, e.g. 4N")
            }
            let n = intbuf.parse::<u32>()?;
            let b = match c.to_ascii_uppercase() {
                CELL_BARCODE_CHAR => ExtractBaseType::CellBarcode,
                UMI_CHAR => ExtractBaseType::Umi,
                MISC_CHAR => ExtractBaseType::Misc,
                e => anyhow::bail!(
                "Illegal character specified for pattern {dstring}: {e}; Choose from C, N, or X."
            ),
            };

            (0..n).for_each(|_| out.push(b.clone()));
            intbuf.clear();
        }
    }

    if !intbuf.is_empty() {
        anyhow::bail!("Trailing integer in pattern string: {intbuf}\nSpecify count before characater, e.g. 4N")
    }

    Ok(out)
}

/// A struct containing information on barcodes extracted from a read given a pattern. This should be cached and overwritten with each record.
#[derive(Clone)]
pub struct Extraction {
    umi_seq: Vec<u8>,
    umi_qual: Vec<u8>,
    cell_seq: Vec<u8>,
    cell_qual: Vec<u8>,
    seq: Vec<u8>,
    seq_qual: Vec<u8>,
    below_quality: bool,
    encoding_delta: u8,
}

impl Extraction {
    fn new(layout: Option<&[ExtractBaseType]>, encoding_delta: u8) -> Self {
        let mut umi_len: usize = 0;
        let mut cell_len: usize = 0;

        if let Some(layout) = layout {
            layout.iter().for_each(|b| match b {
                ExtractBaseType::Umi => umi_len += 1,
                ExtractBaseType::CellBarcode => cell_len += 1,
                _ => (),
            });

            Self {
                umi_seq: Vec::with_capacity(umi_len + 3),
                umi_qual: Vec::with_capacity(umi_len + 3),
                cell_seq: Vec::with_capacity(cell_len + 3),
                cell_qual: Vec::with_capacity(cell_len + 3),
                seq: Vec::new(),
                seq_qual: Vec::new(),
                below_quality: false,
                encoding_delta,
            }
        } else {
            Self {
                umi_seq: Vec::new(),
                umi_qual: Vec::new(),
                cell_seq: Vec::new(),
                cell_qual: Vec::new(),
                seq: Vec::new(),
                seq_qual: Vec::new(),
                below_quality: false,
                encoding_delta,
            }
        }
    }

    fn clear(&mut self) {
        self.umi_seq.clear();
        self.umi_qual.clear();
        self.cell_seq.clear();
        self.cell_qual.clear();
        self.seq.clear();
        self.seq_qual.clear();
        self.below_quality = false;
    }
}

/// Extract a UMI from R1 and R2 (if it exists) given a pattern.
///
/// If a read exists but there is no associated extraction pattern (e.g. R2 exists but pattern2 doesnt), then that read will be unmodified.
///
/// If extraction base phred scores fall below min_qual and mask_below is true, then the base will be replaced with N. Otherwise, it will not be included in the final barcode used to construct the read header.
fn extract_string(
    e: &mut Extraction,
    layout: &[ExtractBaseType],
    read: &Record,
    retain_seq: bool,
    min_qual: u8,
    mask_qual: u8,
) -> Result<(), Error> {
    e.clear();
    let _seq = read.seq();
    let quals = read.qual();
    let mut nbases_extracted = 0;
    let mut qual: u8;

    for (i, b) in layout.iter().enumerate() {
        match b {
            ExtractBaseType::Umi => {
                e.umi_qual.push(quals[i]);

                // dq read
                qual = quals[i] - e.encoding_delta;
                if qual < min_qual {
                    e.below_quality = true;
                    return Ok(());
                // mask crap base
                } else if qual < mask_qual {
                    e.umi_seq.push(b'N');
                // base passes QC
                } else {
                    e.umi_seq.push(_seq[i]);
                }
            }

            ExtractBaseType::CellBarcode => {
                e.cell_seq.push(_seq[i]);
                e.cell_qual.push(quals[i]);
            }

            ExtractBaseType::Misc => {
                e.seq.push(_seq[i]);
                e.seq_qual.push(quals[i]);
            }
        };

        nbases_extracted += 1;
    }

    if retain_seq {
        e.seq.extend(_seq);
        e.seq_qual.extend(read.qual());
    } else {
        e.seq.extend(read.seq()[nbases_extracted..].to_vec());
        e.seq_qual.extend(read.qual()[nbases_extracted..].to_vec());
    }

    Ok(())
}

/// Append a cell and umi barcode, separated by the separator specified to a read header.
///
/// If the read header does not contain a space, it will be appended to the end of the header, otherwise it will be inserted before the last space.
fn remake_read_header(
    cell: &[u8],
    umi: &[u8],
    separator: u8,
    read: &Record,
) -> Result<Vec<u8>, Error> {
    let header = read.id();

    // attempt to append before first space, UMI-tools style, if it exists in the header.
    let (new_qname, _rest) = header.rsplit_once(" ").unwrap_or((header, ""));
    let mut new_qname = new_qname.as_bytes().to_vec();

    if !cell.is_empty() {
        new_qname.push(separator);
        new_qname.extend(cell);
    }

    if !umi.is_empty() {
        new_qname.push(separator);
        new_qname.extend(umi);
    }

    Ok(new_qname)
}

/// Create a new read given header, desc, seq, qual.
fn modify_read(
    header: &[u8],
    desc: Option<&str>,
    seq: &[u8],
    qual: &[u8],
) -> Result<Record, Error> {
    Ok(Record::with_attrs(
        std::str::from_utf8(header)?,
        desc,
        seq,
        qual,
    ))
}

#[derive(Clone)]
/// A struct that stores collections and flags in-between reads processed to avoid constant reallocation.
///
/// It is meant to be cleared after every read, and is expected to be cloned and passed to worker threads working on the same input.
pub struct ExtractionCache {
    layout1: Option<Vec<ExtractBaseType>>,
    layout2: Option<Vec<ExtractBaseType>>,
    e1: Extraction,
    e2: Extraction,
    retain_seq: bool,
    mask_qual: u8,
    min_qual: Option<u8>,
    separator: u8,
}

impl ExtractionCache {
    /// Instantiate from args
    pub fn new_from_args(args: &ExtractArgs) -> Result<Self, Error> {
        let layout1 = if let Some(pattern1) = &args.pattern1 {
            Some(parse_detect_string(pattern1)?)
        } else {
            None
        };

        let layout2 = if let Some(pattern2) = &args.pattern2 {
            Some(parse_detect_string(pattern2)?)
        } else {
            None
        };

        let encoding_delta = match args.qual_encoding {
            FastqQualEncoding::PHRED33 => *PHRED33.start(),
            FastqQualEncoding::PHRED64 => *PHRED64.start(),
            FastqQualEncoding::SOLEXA => *SOLEXA.start(),
        };

        let e1 = Extraction::new(layout1.as_deref(), encoding_delta);
        let e2 = Extraction::new(layout2.as_deref(), encoding_delta);

        Ok(Self {
            layout1,
            layout2,
            e1,
            e2,
            retain_seq: args.retain_seq,
            mask_qual: args.mask_qual,
            min_qual: args.min_qual,
            separator: args.umi_separator as u8,
        })
    }
}

#[derive(Debug)]
struct ThreadScheme {
    workers: usize,
    compress: usize,
}

/// Calculate/distribute available threads to the various methods of the extraction subcommand.
/// We consider:
/// - threads for input
/// - threads for writing AND additional threads for gzip compression
/// - worker threads
///
fn determine_thread_scheme_from_args(args: &ExtractArgs) -> ThreadScheme {
    let paired_input = args.in2.is_some();

    let compress_r1 = args.out1.ends_with(".gz");
    let compress_r2 = args.out2.as_ref().is_some_and(|f| f.ends_with(".gz"));

    let inthreads = if paired_input { 1 } else { 2 };

    let outthreads = if paired_input {
        2 + compress_r1 as i16 + compress_r2 as i16
    } else {
        1 + compress_r1 as i16
    };

    let remainder: i16 = args.threads as i16 - inthreads - outthreads;
    let workers @ compress;

    if remainder <= 0 {
        workers = 1;
        compress = 1;
    } else {
        workers = (remainder as f32 * 0.6).ceil() as usize;
        compress = (remainder as f32 * 0.4).ceil() as usize;
    }

    ThreadScheme { workers, compress }
}

/// The main extraction function that operates on a read pair, extracting barcodes, concatenating them in the case where R2 and P2 are supplied, and adding them to the read header.
fn process_pair2(mut pair: ReadPair, ecache: &mut ExtractionCache) -> Result<ReadPair, Error> {
    // we'll be merging R1 and R2 barcodes if both exist, so double cap
    let mut umi: Vec<u8> = Vec::with_capacity(ecache.e1.umi_seq.capacity() * 2);
    let mut cell: Vec<u8> = Vec::with_capacity(ecache.e1.cell_seq.capacity() * 2);

    match (&ecache.layout1, pair.r1.as_ref()) {
        // we either have no pattern or no r1, so do nothing
        (None, _) | (_, None) => (),

        // process
        (Some(pattern1), Some(r1)) => {
            extract_string(
                &mut ecache.e1,
                pattern1,
                r1,
                ecache.retain_seq,
                ecache.min_qual.unwrap_or(0),
                ecache.mask_qual,
            )?;

            umi.extend(&ecache.e1.umi_seq);
            cell.extend(&ecache.e1.cell_seq);
        }
    }

    if ecache.e1.below_quality {
        pair.r1 = None;
        pair.r2 = None;
        return Ok(pair);
    };

    match (&ecache.layout2, pair.r2.as_ref()) {
        // we either have no pattern or no r2, so do nothing
        (None, _) | (_, None) => (),

        // process
        (Some(pattern2), Some(r2)) => {
            extract_string(
                &mut ecache.e2,
                pattern2,
                r2,
                ecache.retain_seq,
                ecache.min_qual.unwrap_or(0),
                ecache.mask_qual,
            )?;

            umi.extend(&ecache.e2.umi_seq);
            cell.extend(&ecache.e2.cell_seq);
        }
    }

    if ecache.e2.below_quality {
        pair.r1 = None;
        pair.r2 = None;
        return Ok(pair);
    };

    // Re-header R1 if it exists
    if let Some(r1) = pair.r1.as_mut() {
        if !ecache.e1.below_quality {
            let header = remake_read_header(&cell, &umi, ecache.separator, r1)?;
            pair.r1 = Some(modify_read(
                &header,
                // r1.desc(),
                None,
                &ecache.e1.seq,
                &ecache.e1.seq_qual,
            )?);
        } else {
            pair.r1 = None;
        }
    }

    // Re-header R2 if it exists
    if let Some(r2) = pair.r2.as_mut() {
        if !ecache.e2.below_quality {
            let header = remake_read_header(&cell, &umi, ecache.separator, r2)?;
            pair.r2 = Some(modify_read(
                &header,
                // r2.desc(),
                None,
                &ecache.e2.seq,
                &ecache.e2.seq_qual,
            )?);
        } else {
            pair.r2 = None;
        }
    }

    Ok(pair)
}

#[derive(Debug)]

/// A single worker to extract barcodes from reads.
/// It maintains a handle to its channel  for recieving reads, and a joinhandle.
pub struct ExtractionWorker {
    input_handle: Option<Sender<IntakeOrdered>>,
    handle: Option<std::thread::JoinHandle<Result<(), Error>>>,
}

impl ExtractionWorker {
    /// Spin up the extraction worker as its own thread, and process reads based off the [ExtractionCache] supplied.
    ///
    /// A `Sender` must be supplied; the thread will send processed reads to this sender for downstream operations.
    fn run(&mut self, write_handle: Sender<IntakeOrdered>, mut ecache: ExtractionCache) {
        let (s, r): (Sender<IntakeOrdered>, Receiver<IntakeOrdered>) = unbounded();

        let handle = std::thread::spawn(move || -> Result<(), Error> {
            let mut reads_processed: usize = 0;
            while let Ok(mut rec) = r.recv() {
                let r = process_pair2(rec.pair, &mut ecache)?;
                rec.pair = r;
                write_handle.send(rec)?;
                reads_processed += 1;
            }

            println! {"Reads processed: {reads_processed}"};

            Ok(())
        });

        self.handle = Some(handle);
        self.input_handle = Some(s);
    }

    /// Instantiate an extraction worker, but do not spin up a thread or link to an output.
    fn new() -> Self {
        Self {
            input_handle: None,
            handle: None,
        }
    }
}

/// A simple struct that counts number of reads input, maintains a list of extraction workers.
pub struct ExtractionWorkerPool {
    workers: Vec<ExtractionWorker>,
    n_workers: usize,
    n_loaded: usize,
}

impl ExtractionWorkerPool {
    /// Create a thread pool with a given number of workers.
    fn new(n_workers: usize) -> Self {
        assert!(n_workers > 0);

        let mut workers = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            workers.push(ExtractionWorker::new());
        }

        Self {
            workers,
            n_workers,
            n_loaded: 0,
        }
    }

    /// Spin up each worker thread to idle and wait for input reads to process and send.
    ///
    /// An output Sender must be supplied, so each worker can offload processed reads.
    fn start(&mut self, write_handle: Sender<IntakeOrdered>, ecache: &ExtractionCache) {
        for worker in self.workers.iter_mut() {
            worker.run(write_handle.clone(), ecache.clone());
        }
    }

    /// Accept a read into the thread pool, sending it to one of the worker threads and recording its order in the file being read.
    ///
    /// Currently, reads are distributed across workers in a round-robin fashion.
    fn intake_read(&mut self, r: ReadPair) {
        let idx = self.n_loaded % self.n_workers;

        let intake = IntakeOrdered {
            pair: r,
            order: self.n_loaded as u32,
        };

        self.workers[idx]
            .input_handle
            .as_ref()
            .unwrap()
            .send(intake)
            .unwrap();

        self.n_loaded += 1;
    }

    /// Signal all worker threads to finish by dropping their input handles and wait for them finish.
    fn end(&mut self) -> Result<(), Error> {
        for worker in self.workers.drain(..) {
            drop(worker.input_handle);

            worker
                .handle
                .context("Nonexisting worker handle")?
                .join()
                .expect("Failed to join worker handle")?;
        }

        Ok(())
    }
}

pub fn run_extract(args: &ExtractArgs) -> Result<(), Error> {
    let threads = determine_thread_scheme_from_args(args);

    let mut fio = FastqIO::new_from_inputs(
        &args.in1,
        &args.in2,
        &args.out1,
        &args.out2,
        args.batch_size,
    )?;
    fio.start(threads.compress)?;

    let ecache = ExtractionCache::new_from_args(args)?;

    let mut pool = ExtractionWorkerPool::new(threads.workers);
    pool.start(fio.get_output_handle()?, &ecache);

    loop {
        match fio.next_pair() {
            Err(e) => anyhow::bail!("Error reading fastq input: {e}"),
            Ok(next) => match next {
                Some(pair) => pool.intake_read(pair),
                None => break, // end of file
            },
        }
    }

    pool.end()?;
    fio.terminate()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_numerical_detect_string1() {
        let dstring = "4N7C";
        let out = parse_detect_string_alphanumeric(dstring);

        let expected = std::iter::repeat_n(ExtractBaseType::Umi, 4)
            .chain(std::iter::repeat_n(ExtractBaseType::CellBarcode, 7))
            .collect::<Vec<ExtractBaseType>>();

        assert!(out.is_ok());
        assert_eq!(out.unwrap(), expected)
    }

    #[test]
    fn test_numerical_detect_string2() {
        let dstring = "N47C";
        let out = parse_detect_string_alphanumeric(dstring);

        assert!(out.is_err());
    }

    #[test]
    fn test_numerical_detect_string3() {
        let dstring = "4N7C8";
        let out = parse_detect_string_alphanumeric(dstring);

        assert!(out.is_err());
    }
}
