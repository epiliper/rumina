use crate::cli::{dedup_args::DedupArgs, extract_args::ExtractArgs};
use clap::{Parser, Subcommand};

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Deduplicate or cluster reads based on UMI barcodes, with error correction.
    Dedup(DedupArgs),
    /// Extract UMI barcodes from read sequence in fastq/fastq.gz files.
    Extract(ExtractArgs),

    Test,
}

#[derive(Debug, Parser)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,
}
