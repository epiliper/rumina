use crate::cli::DedupArgs;
use anyhow::Error;

pub trait FileProcess {
    fn init_from_args(
        args: &DedupArgs,
        bam_file_path: &str,
        bam_file_name: &str,
    ) -> Result<Self, Error>
    where
        Self: Sized;

    fn process(self) -> Result<(), Error>;
}
