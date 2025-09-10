// use crate::fastq_output::FastqOutput;
// use crate::fastqio::{FastqInput, IntakeOrdered, ReadPair};
use crate::io::{
    fastq_output::FastqOutput,
    fastqio::{FastqInput, IntakeOrdered, ReadPair},
};
use anyhow::{Context, Error};
use crossbeam::channel::Sender;
use std::path::{Path, PathBuf};

/// An abstraction over readers and writers to either one or two fastq inputs, expected to be R1 and R2.
/// See `FastqInput` and `FastqOutput` for more details.
pub struct FastqIO {
    r1: FastqInput,
    r2: Option<FastqInput>,
    out: FastqOutput,
    pub writer_out_handle: Option<Sender<IntakeOrdered>>,
}

impl FastqIO {
    pub fn new_from_inputs(
        in1: &Path,
        in2: &Option<PathBuf>,
        out1: &Path,
        out2: &Option<PathBuf>,
        cache_size: usize,
    ) -> Result<Self, Error> {
        let r1 = FastqInput::from_file(in1)?;

        let r2: Option<FastqInput> = if let Some(in2) = in2 {
            if out2.is_none() {
                anyhow::bail!("Must supply output for file R2!");
            }

            Some(FastqInput::from_file(in2)?)
        } else {
            None
        };

        let out = FastqOutput::init_from_outputs(out1, out2, cache_size)?;

        Ok(Self {
            r1,
            r2,
            out,
            writer_out_handle: None,
        })
    }

    pub fn get_output_handle(&self) -> Result<Sender<IntakeOrdered>, Error> {
        self.writer_out_handle
            .as_ref()
            .context("Writer channel unitialized!")
            .cloned()
    }

    /// Read from R1 and R2 (if it was supplied), consoliding read errors from either into one error.
    pub fn next_pair(&mut self) -> Result<Option<ReadPair>, Error> {
        let ret1 = self
            .r1
            .records
            .next_record()
            .transpose()
            .context("Error reading R1!")?;

        let ret2;

        if let Some(r2) = self.r2.as_mut() {
            ret2 = r2
                .records
                .next_record()
                .transpose()
                .context("Error reading R2!")?;
        } else {
            ret2 = None;
        }

        if ret1.is_none() && ret2.is_none() {
            Ok(None)
        } else {
            Ok(Some(ReadPair { r1: ret1, r2: ret2 }))
        }
    }

    pub fn start(&mut self, out_threads: usize) -> Result<(), Error> {
        self.out.run(out_threads)?;
        self.writer_out_handle = self.out.input_handle.as_ref().cloned();
        Ok(())
    }

    pub fn terminate(self) -> Result<(), Error> {
        drop(self.writer_out_handle);
        self.out.terminate()
    }
}
