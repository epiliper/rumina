use crate::args::Args;
use anyhow::Error;

pub trait FileProcess {
    fn init_from_args(
        args: &Args,
        bam_file_path: &String,
        bam_file_name: &String,
    ) -> Result<Self, Error>
    where
        Self: Sized;

    fn process(self) -> Result<(), Error>;
}
