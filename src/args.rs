use clap::{Parser, ValueEnum};
use colored::Colorize;

#[derive(ValueEnum, Debug, Clone)]
pub enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}
#[derive(Parser, Debug)]
#[command(version, about, term_width = 0)]
pub struct Args {
    #[arg(short = 'i', index = 1)]
    /// a BAM file or directory of BAM files to be processed.
    /// BAM files must have UMIs in the read QNAME, and be sorted and indexed.
    pub input: String,

    #[arg(short = 'g', long = "grouping_method")]
    /// specifies UMI merging method for UMI error correction
    /// acyclic: Same as directional, but networks are limited to depth of one
    /// raw: Treat each raw UMI as genuine; no merging
    pub grouping_method: GroupingMethod,

    #[arg(short = 's', long = "separator")]
    /// character in read QNAME delimiting UMI barcode. Usually '_' or '.'
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
    /// Use only R1 for deduplication, discard R1, similar to UMI-tools
    pub r1_only: bool,

    #[arg(short = 't', long = "threads", default_value_t = num_cpus::get())]
    /// number of threads to use. Will default to use all if not specified.
    pub threads: usize,

    #[arg(long = "strict_threads", default_value_t = false)]
    /// Also constrain SAM read/write operations to use the number of threads specified in
    /// --threads. Choose this option to make CPU usage more predictable, but expect runtime to
    /// increase when using fewer than all available threads.
    pub strict_threads: bool,
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
            {}: {}\n",
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
            "Singletons".purple(),
            self.singletons,
            "Halve pairs".purple(),
            self.r1_only
        )?;

        Ok(())
    }
}

pub fn parse_args() -> Args {
    Args::parse()
}
