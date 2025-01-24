use clap::{Parser, Subcommand, ValueEnum};

#[derive(ValueEnum, Debug, Clone)]
pub enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}

#[derive(Parser, Debug)]
#[command(term_width = 0)]
pub struct Args {
    #[command(subcommand)]
    pub command: Command,

    #[arg(long = "in", index = 1)]
    pub input: String,

    #[arg(long = "out", index = 2)]
    pub output: String,

    #[arg(long = "threads", index = 3)]
    pub threads: usize,

    #[arg(long = "split_window")]
    pub split_window: Option<i64>,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    Merge {
        #[arg(long = "merge_pairs")]
        ref_fasta: String,

        #[arg(long = "min_overlap_bp")]
        min_overlap_bp: usize,
    },

    Group {
        #[arg(short = 's')]
        separator: String,

        #[arg(short = 'g')]
        grouping_method: GroupingMethod,

        #[arg(long = "length")]
        length: bool,

        #[arg(long = "only-group")]
        only_group: bool,

        #[arg(long = "singletons")]
        singletons: bool,

        #[arg(long = "r1-only")]
        r1_only: bool,
    },
}

pub fn parse_args() -> Args {
    Args::parse()
}
