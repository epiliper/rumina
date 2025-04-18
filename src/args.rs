use clap::{Parser, ValueEnum};
use colored::Colorize;

#[derive(ValueEnum, Debug, Clone)]
pub enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}
#[derive(Parser, Debug)]
#[command(version, about)]
pub struct Args {
    #[arg(index = 1)]
    /// a BAM file or directory of BAM files to be processed.
    /// BAM files must have UMIs in the read QNAME, and be sorted and indexed.
    pub input: String,

    #[arg(short = 'g', long = "grouping_method")]
    /// specifies UMI merging method for UMI error correction
    /// directional (standard): predict mutant UMIs based on hamming distance from other UMIs
    /// acyclic: Same as directional, but networks are limited to depth of one
    /// raw: Treat each raw UMI as genuine; no merging
    pub grouping_method: GroupingMethod,

    #[arg(short = 's', long = "separator")]
    /// character in read QNAME delimiting UMI barcode. Usually '_' or ':'
    pub separator: String,

    #[arg(long = "outdir", default_value = "rumina_output")]
    /// directory (relative to parent dir of input file) in which to store output files. Will be created if it doesn't exist
    pub outdir: String,

    #[arg(value_parser = clap::value_parser!(i64).range(1..), short = 'x', long = "split_window")]
    /// BAM file splitting strategy. Not using this arg or passing in 0 will have all BAM
    /// coordinates processed at once, which may be faster but also incur significant memory usage
    /// with larger alignments.
    pub split_window: Option<i64>,

    #[arg(long = "merge_pairs", conflicts_with = "halve_pairs")]
    /// merge overlapping forward/reverse read pairs with the same UMI. Requires reference genome
    /// FASTA for realignment.
    pub merge_pairs: Option<String>,

    #[arg(long = "min_overlap_bp", default_value = "3", value_parser = clap::value_parser!(i64).range(1..))]
    /// minimum bases for read overlap. See "--merge_pairs"
    pub min_overlap_bp: i64,

    #[arg(short = 'l', long = "length")]
    /// group reads by length in addition to reference coordinate
    pub length: bool,

    #[arg(long = "only-group")]
    /// only group reads; do not deduplicate
    pub only_group: bool,

    #[arg(long = "singletons")]
    /// if used, singleton groups (1-2 read UMI clusters) will be processed like other groups,
    /// instead of discarded.
    pub singletons: bool,

    #[arg(long = "halve_pairs", conflicts_with = "merge_pairs")]
    /// Use only R1 for deduplication, discard R2, similar to UMI-tools
    pub r1_only: bool,

    #[arg(short = 't', long = "threads", default_value_t = num_cpus::get())]
    /// number of threads to use. Will default to use all if not specified.
    pub threads: usize,

    #[arg(long = "strict_threads", default_value_t = false)]
    /// Also constrain SAM read/write operations to use the number of threads specified in
    /// --threads. Choose this option to make CPU usage more predictable, but expect runtime to
    /// increase when using fewer than all available threads.
    pub strict_threads: bool,

    #[arg(short = 'p', long = "percentage", default_value_t = 0.5)]
    /// The percentage of a parent UMIs read count an offshoot must hold, at maximum, to be
    /// considered an offshoot. This is one of two criteria for clustering; also see -m /
    /// --max_edit.
    pub percentage: f32,

    #[arg(short = 'm', long = "max_edit", default_value_t = 1, value_parser = clap::value_parser!(u32).range(0..))]
    /// The maximum hamming distance difference between two UMIs for one to be considered an
    /// offshoot of the other.
    pub max_edit: u32,
}

impl std::fmt::Display for Args {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}: {}\n\
            {}: {:?}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {:?}\n\
            {}: {:?}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n
",
            "Input".purple(),
            self.input,
            "Grouping method".purple(),
            self.grouping_method,
            "Separator".purple(),
            self.separator,
            "Outdir".purple(),
            self.outdir,
            "Threads".purple(),
            self.threads,
            "Split window".purple(),
            self.split_window,
            "Merge pairs".purple(),
            self.merge_pairs,
            "Min overlap BP".purple(),
            self.min_overlap_bp,
            "Length".purple(),
            self.length,
            "Group only".purple(),
            self.only_group,
            "Keep singletons".purple(),
            self.singletons,
            "Halve pairs".purple(),
            self.r1_only,
            "Percentage".purple(),
            self.percentage,
            "Max edit distance".purple(),
            self.max_edit,
        )?;

        Ok(())
    }
}

pub fn parse_args() -> Args {
    let args = Args::parse();
    if args.percentage < 0.01 || args.percentage > 1.0 {
        panic!(
            "Invalid value {} for -p/--percentage! Choose a value between (inclusive) 0.01 and 1.0",
            args.percentage
        );
    }
    args
}
