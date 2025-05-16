use anyhow::Error;
use clap::{Parser, ValueEnum};
use colored::Colorize;
use indoc::{formatdoc, indoc};

const DEFAULT_PERCENT: f32 = 0.5;
const DEFAULT_MIN_DEPTH: usize = 3;
const DEFAULT_MAX_EDIT: u32 = 1;
const DEFAULT_MIN_OVERLAP: i64 = 3;
static DEFAULT_THREADS: std::sync::LazyLock<usize> = std::sync::LazyLock::new(|| num_cpus::get());

#[derive(ValueEnum, Debug, Clone)]
pub enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}
#[derive(Parser, Debug)]
#[command(version, about)]
pub struct Args {
    #[arg(index = 1, help_heading = "REQUIRED ARGUMENTS")]
    /// a BAM file or directory of BAM files to be processed.
    /// BAM files must have UMIs in the read QNAME, and be sorted and indexed.
    pub input: String,

    #[arg(
        short = 'g',
        long = "grouping_method",
        hide_possible_values = true,
        help_heading = "REQUIRED ARGUMENTS",
        help = indoc! {"
            Specifes UMI clustering method: choices:
                - directional: as in UMI-tools, predict mutant UMIs based on hamming distance and frequency from other UMIs
                - acyclic: same as directional, but networks are limited to depth of one; recommended for low-PCR datasets with longer barcodes
                - raw: treat UMIs as-is; no grouping. Use with caution.
            "}
    )]
    pub grouping_method: GroupingMethod,

    #[arg(short = 's', long = "separator", help_heading = "REQUIRED ARGUMENTS", help = indoc! {"
        The character in the read QNAME delimiting the UMI barcode from the rest of the QNAME. Usually '_' or ':'
        "})]
    pub separator: String,

    #[arg(short = 'p', long = "percentage", default_value_t = DEFAULT_PERCENT, help_heading = "CLUSTERING OPTIONS", hide_default_value = true, help = formatdoc!{"
        The fraction of a parent UMI's read count an offshoot's count must be, at maximum, to be considered an offshoot [{default}]
        ", default = DEFAULT_PERCENT})]
    pub percentage: f32,

    #[arg(short = 'm', long = "max_edit", default_value_t = DEFAULT_MAX_EDIT, hide_default_value = true, value_parser = clap::value_parser!(u32).range(0..), help_heading = "CLUSTERING OPTIONS", help = formatdoc!{"
        The maximum hamming distance difference between two UMIs for one to be considered an immediate offshoot of the other [{default}]
        ", default = DEFAULT_MAX_EDIT})]
    pub max_edit: u32,

    #[arg(short = 'd', long = "min_depth", default_value_t = DEFAULT_MIN_DEPTH, hide_default_value = true, help_heading = "CLUSTERING OPTIONS", help = formatdoc! {
        "Minimum number of reads in a cluster for it be processed and output, overwritten by --singletons [{default}]
        ", default = DEFAULT_MIN_DEPTH})]
    pub min_cluster_depth: usize,

    #[arg(short = 'u', long = "rev", help_heading = "CLUSTERING OPTIONS", help = indoc! {"
        Search also for reverse complement barcode when clustering"})]
    pub cluster_rev: bool,

    #[arg(short = 'l', long = "length", help_heading = "CLUSTERING OPTIONS", help = indoc!{"group reads by length in addition to reference coordinate before clustering.
        "})]
    pub length: bool,

    #[arg(short = 'v', long = "only-group", help_heading = "CLUSTERING OPTIONS", help = indoc!{" only group reads; do not deduplicate; useful for manually inspecting grouping."})]
    pub only_group: bool,

    #[arg(short = 'f', long = "singletons", help_heading = "CLUSTERING OPTIONS", help = indoc!{"
        remove minimum depth limit for UMI clusters to be processed; identical to --min_depth 1
        "})]
    pub singletons: bool,

    #[arg(short = 'o', long = "outdir", default_value = "rumina_output", help_heading = "MISC OPTIONS", help = indoc! {"
        directory (relative to parent dir of input) in which to store output files
        "})]
    pub outdir: String,

    #[arg(value_parser = clap::value_parser!(i64).range(1..), short = 'x', long = "split_window", help_heading = "PERFORMANCE OPTIONS", help = indoc!{"
        number of reference coordinates of a given BAM to process at a time [no limit]; a starting value of 100 is recommended if memory is limited"})]
    pub split_window: Option<i64>,

    #[arg(short = 'm', long = "merge_pairs", conflicts_with = "halve_pairs", help_heading = "PAIRED-END OPTIONS", help = indoc!{"
        Use the supplied reference fasta to merge overlapping forward/reverse read pairs with the same UMI. Only supports single-fasta files currently. For an alternative paired-end strategy, see --halve_pairs.
        "})]
    pub merge_pairs: Option<String>,

    #[arg(short = 'b', long = "min_overlap_bp", default_value_t = DEFAULT_MIN_OVERLAP, hide_default_value = true, value_parser = clap::value_parser!(i64).range(1..), help_heading = "PAIRED-END OPTIONS", help = formatdoc!{"
        The minimum number of overlapping bases for two reads to be merged when using --merge_pairs [{default}]
        ", default = DEFAULT_MIN_DEPTH})]
    pub min_overlap_bp: i64,

    #[arg(short = 'l', long = "halve_pairs", conflicts_with = "merge_pairs", help_heading = "PAIRED-END OPTIONS", help = indoc!{"
        Use only R1 for deduplication, and pair R1 with deduplicated R2 for output, similar to UMI-tools
        "})]
    pub r1_only: bool,

    #[arg(short = 't', long = "threads", default_value_t = *DEFAULT_THREADS, hide_default_value = true, help_heading = "PERFORMANCE OPTIONS", help = formatdoc! {"
        Number of threads to use for parallel processing of coordinates. Will default to all if not specified [{}].
        ", *DEFAULT_THREADS})]
    pub threads: usize,

    #[arg(short = 'c', long = "strict_threads", default_value_t = false, help_heading = "PERFORMANCE OPTIONS", help = indoc! {"Constrain sam read/write operations to --threads
        "})]
    pub strict_threads: bool,

    #[arg(short = 'q', long = "progress", default_value_t = false, help_heading = "MISC OPTIONS", help = indoc!{"Show progress bars for read clustering information."})]
    pub progress: bool,
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

pub fn parse_args() -> Result<Args, Error> {
    let args = Args::parse();
    if args.percentage < 0.01 || args.percentage > 1.0 {
        anyhow::bail!(
            "Invalid value {} for -p/--percentage! Choose a value between (inclusive) 0.01 and 1.0",
            args.percentage
        )
    }
    Ok(args)
}
