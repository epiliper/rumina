use crate::io::fastqio::{
    fastq_create_writer_compressed, fastq_create_writer_decompressed, WritesFastqRecords,
};
use anyhow::{Context, Error};
use bio::io::fastq::Record;
use crossbeam::channel::{unbounded, Receiver, Sender};
use std::boxed::Box;
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

/// A struct that synchronizes input reads from multiple threads into a single collection fed into a Fastq writer for output to a file or stdout.
pub struct ChanneledFastqWriter {
    join_handle: Option<JoinHandle<Result<u32, Error>>>,
    inner: InnerFastqWriter,
    pub input_handle: Option<Sender<Record>>,
}

/// A struct that launches fastq writer(s) on their own threads, each with their own input channels.
struct InnerFastqWriter {
    chunk_size: usize, // how many reads to accumulate before writing
    outfile: PathBuf,  // the output file
}

pub fn create_writer_from_file(
    infile: &Path,
    threads: usize,
) -> Result<Box<dyn WritesFastqRecords>, Error> {
    let is_gzip = match infile.to_str().unwrap().rsplit_once(".") {
        None => anyhow::bail!("Unable to detect file format!"),
        Some((_pre, end)) => end == "gz",
    };

    let writer: Box<dyn WritesFastqRecords>;
    if is_gzip {
        writer = Box::new(fastq_create_writer_compressed(infile, threads)?);
    } else {
        writer = Box::new(fastq_create_writer_decompressed(infile, threads)?);
    }

    Ok(writer)
}

impl InnerFastqWriter {
    /// Given an output file, spawn a thread that accumulates input reads via a channel and writes them to output once they accumulate to an arbitrary quantity.
    #[allow(clippy::complexity)]
    fn run(
        &self,
        threads: usize,
    ) -> Result<(Sender<Record>, JoinHandle<Result<u32, Error>>), Error> {
        let mut buffer: Vec<Record> = Vec::with_capacity(self.chunk_size);
        let chunk_size = self.chunk_size;
        let outfile = self.outfile.clone();

        let is_gzip = match self.outfile.to_str().unwrap().rsplit_once(".") {
            None => anyhow::bail!("Unable to detect file format!"),
            Some((_pre, end)) => end == "gz",
        };

        let (s, r): (Sender<Record>, Receiver<Record>) = unbounded();

        Ok((
            s,
            std::thread::spawn(move || -> Result<u32, Error> {
                let mut writer: Box<dyn WritesFastqRecords>;
                if is_gzip {
                    writer = Box::new(fastq_create_writer_compressed(outfile.as_path(), threads)?);
                } else {
                    writer = Box::new(fastq_create_writer_decompressed(
                        outfile.as_path(),
                        threads,
                    )?);
                }
                let mut num_written = 0;
                while let Ok(r) = r.recv() {
                    buffer.push(r);
                    if buffer.len() >= chunk_size {
                        buffer.drain(..).for_each(|r| {
                            writer.write_record(r).expect("bad read write");
                            num_written += 1;
                        })
                    }
                }

                // Sender dropped. Time to wrap-up and go home.
                buffer
                    .drain(..)
                    .for_each(|r| writer.write_record(r).expect("bad read write"));

                Ok(num_written)
            }),
        ))
    }
}

impl ChanneledFastqWriter {
    pub fn create_from_outfile(outfile: &Path, cap: Option<usize>) -> Result<Self, Error> {
        let inner = InnerFastqWriter {
            chunk_size: cap.unwrap_or(1e4 as usize),
            outfile: outfile.to_path_buf(),
        };

        Ok(Self {
            inner,
            join_handle: None,
            input_handle: None,
        })
    }

    /// Spin up the writer thread.
    pub fn run(&mut self, threads: usize) -> Result<(), Error> {
        let (input_handle, join_handle) = self.inner.run(threads)?;
        self.input_handle = Some(input_handle);
        self.join_handle = Some(join_handle);

        Ok(())
    }

    /// Signal for the writer thread to end, and wait for it to do so.
    pub fn end(mut self) -> Result<(), Error> {
        drop(self.input_handle);
        self.join_handle
            .take()
            .context("Writer thread unitialized")?
            .join()
            .expect("Failed to join writer thread")?;

        Ok(())
    }
}
