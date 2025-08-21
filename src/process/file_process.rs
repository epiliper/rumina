use crate::args::Args;
use anyhow::Error;

pub trait FileProcess {
    fn init_from_args(args: &Args, bam_file_path: &str, bam_file_name: &str) -> Result<Self, Error>
    where
        Self: Sized;

    fn process(self) -> Result<(), Error>;
}
