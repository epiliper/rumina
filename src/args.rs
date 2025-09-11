use crate::cli::{dedup_args::DedupArgs, extract_args::ExtractArgs};
use anyhow::Error;
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Subcommand)]
pub enum Command {
    Dedup(DedupArgs),
    Group(DedupArgs),
    Extract(ExtractArgs),
}

#[derive(Debug, Parser)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,
}

pub fn parse_args_dedup() -> Result<DedupArgs, Error> {
    let args = DedupArgs::parse();
    if args.percentage < 0.01 || args.percentage > 1.0 {
        anyhow::bail!(
            "Invalid value {} for -p/--percentage! Choose a value between (inclusive) 0.01 and 1.0",
            args.percentage
        )
    }
    Ok(args)
}

pub fn parse_args_extract() -> Result<ExtractArgs, Error> {
    let args = ExtractArgs::parse();
    Ok(args)
}
