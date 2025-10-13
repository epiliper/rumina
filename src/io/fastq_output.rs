use crate::io::fastq_writer::ChanneledFastqWriter;

use anyhow::{Context, Error};
use crossbeam::channel::{unbounded, Receiver, Sender};
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

use crate::io::fastqio::IntakeOrdered;

static MIN_WRITE_THRESHOLD: usize = 100;

/// A struct that maintains writers R1, and R2 if it is needed. It consolidates input to both writers into a single channel, where each channel item is a `ReadPair`.
///
/// Because inputs arrive as R1/R2 pairs and are processed sequentially, the pairedness of outputs is maintained.  
pub struct FastqOutput {
    r1: ChanneledFastqWriter,
    r2: Option<ChanneledFastqWriter>,
    join_handle: Option<JoinHandle<Result<(), Error>>>,
    pub input_handle: Option<Sender<IntakeOrdered>>,
    batch_size: usize,
}

impl FastqOutput {
    /// Create writers for supplied outputs, but don't spin up their threads yet.
    pub fn init_from_outputs(
        out1: &Path,
        out2: &Option<PathBuf>,
        batch_size: usize,
    ) -> Result<Self, Error> {
        if !Path::exists(Path::parent(out1).context("Unable to resolve path for R1 output")?) {
            anyhow::bail!(
                "At least one level of output directory for R1 does not exist in file system: {:?}",
                out1
            );
        }
        let r1 = ChanneledFastqWriter::create_from_outfile(out1, Some(batch_size))?;

        let r2 = if let Some(out2) = out2 {
            if !Path::exists(Path::parent(out2).context("Unable to resolve path for R2 output")?) {
                anyhow::bail!(
                    "At least one level of output directory for R2 does not exist in file system: {:?}",
                    out1
                );
            }
            Some(ChanneledFastqWriter::create_from_outfile(
                out2,
                Some(batch_size),
            )?)
        } else {
            None
        };

        Ok(Self {
            r1,
            r2,
            join_handle: None,
            input_handle: None,
            batch_size,
        })
    }

    /// Spin up the threads for all writers.
    ///
    /// Items in the output channel to be written are sorted by their read-in order, and writing to file only proceeds while items are sequential, e.g.:
    ///
    /// pair1.order = 4, pair2.order = 5, pair3.order = 6
    ///
    /// We stop writing once i.order > (i - 1).order + 1, and wait for more inputs.
    ///
    /// Order here is intended to be the order in which reads were read in from an input file, but can be any number.
    pub fn run(&mut self, threads: usize) -> Result<(), Error> {
        let (s, r): (Sender<IntakeOrdered>, Receiver<IntakeOrdered>) = unbounded();

        let r1_threads = std::cmp::min(threads / (1 + self.r2.is_some() as usize), 1);
        let r2_threads = self.r2.is_some() as usize;

        // spin up r1 writer
        self.r1.run(r1_threads)?;

        let in1 = self
            .r1
            .input_handle
            .take()
            .context("R1 writer not initialized!")?;

        // spin up r2 writer
        let mut in2 = if let Some(r2) = self.r2.as_mut() {
            r2.run(r2_threads)?;
            Some(
                r2.input_handle
                    .take()
                    .context("R2 writer not initialized!")?,
            )
        } else {
            None
        };

        let batch_size = self.batch_size;

        let write_threshold =
            std::cmp::min((self.batch_size as f32 * 0.9) as usize, MIN_WRITE_THRESHOLD);

        self.join_handle = Some(std::thread::spawn(move || -> Result<(), Error> {
            let mut buffer: Vec<IntakeOrdered> = Vec::with_capacity(batch_size);
            let mut next_expected_order = 0;

            while let Ok(r) = r.recv() {
                buffer.push(r);

                if buffer.len() >= write_threshold {
                    buffer.sort_by(|a, b| a.order.cmp(&b.order));

                    let mut processable_count = 0;
                    for item in &buffer {
                        if item.order == next_expected_order {
                            processable_count += 1;
                            next_expected_order += 1;
                        } else {
                            // out of order, so stop writing and wait for more inputs
                            break;
                        }
                    }

                    for i in buffer.drain(..processable_count) {
                        if let Some(r1) = i.pair.r1 {
                            in1.send(r1)?;
                        }

                        if let Some(r2) = i.pair.r2 {
                            in2.as_mut().context("no R2 reader")?.send(r2)?;
                        }
                    }
                }
            }

            buffer.sort_by(|a, b| a.order.cmp(&b.order));

            for i in buffer.drain(..) {
                if let Some(r1) = i.pair.r1 {
                    in1.send(r1)?;
                }

                if let Some(r2) = i.pair.r2 {
                    in2.as_mut().context("no R2 reader")?.send(r2)?;
                }
            }

            Ok(())
        }));

        self.input_handle = Some(s);

        Ok(())
    }

    pub fn terminate(mut self) -> Result<(), Error> {
        // terminate read pair channel
        drop(self.input_handle);

        // wait for thread to finish sending to writers
        self.join_handle
            .take()
            .context("Writer uninitialized")?
            .join()
            .expect("Failed to join writer thread")?;

        // terminate the writers themselves
        self.r1.end()?;

        if let Some(r2) = self.r2 {
            r2.end()?;
        }

        Ok(())
    }
}
