use anyhow::Error;
use clap::{Parser, ValueEnum};
use colored::Colorize;

const DEFAULT_PERCENT: f32 = 0.5;
const DEFAULT_MIN_DEPTH: usize = 3;
const DEFAULT_MAX_EDIT: u32 = 1;
const DEFAULT_MIN_OVERLAP: i64 = 3;
static DEFAULT_THREADS: std::sync::LazyLock<usize> = std::sync::LazyLock::new(num_cpus::get);

#[derive(ValueEnum, Debug, Clone)]
pub enum GroupingMethod {
    Acyclic,
    Directional,
    Raw,
}
#[derive(Parser, Debug)]
#[command(version, about, override_help = DEDUP_HELP)]
pub struct DedupArgs {
    #[clap(short = 'i')]
    pub input: String,

    #[arg(short = 'g', long = "grouping_method")]
    pub grouping_method: GroupingMethod,

    #[arg(short = 's', long = "separator")]
    pub separator: String,

    #[arg(short = 'p', long = "percentage", default_value_t = DEFAULT_PERCENT)]
    pub percentage: f32,

    #[arg(short = 'm', long = "max-edit", default_value_t = DEFAULT_MAX_EDIT)]
    pub max_edit: u32,

    #[arg(short = 'd', long = "min-depth", default_value_t = DEFAULT_MIN_DEPTH)]
    pub min_cluster_depth: usize,

    #[arg(short = 'u', long = "rev")]
    pub cluster_rev: bool,

    #[arg(short = 'l', long = "length")]
    pub length: bool,

    #[arg(short = 'v', long = "only-group")]
    pub only_group: bool,

    #[arg(short = 'f', long = "singletons")]
    pub singletons: bool,

    #[arg(short = 'o', long = "outdir", default_value = "rumina_output")]
    pub outdir: String,

    #[arg(long = "sort", help_heading = "MISC OPTIONS")]
    pub ensure_sorted: bool,

    #[arg(value_parser = clap::value_parser!(i64).range(1..), short = 'x', long = "split-window")]
    pub split_window: Option<i64>,

    #[arg(short = 'm', long = "merge-pairs", conflicts_with = "paired")]
    pub merge_pairs: Option<String>,

    #[arg(short = 'b', long = "min-overlap-bp", default_value_t = DEFAULT_MIN_OVERLAP)]
    pub min_overlap_bp: i64,

    #[arg(short = 'l', long = "paired", conflicts_with = "merge_pairs")]
    pub paired: bool,

    #[arg(short = 't', long = "threads", default_value_t = *DEFAULT_THREADS)]
    pub threads: usize,

    #[arg(short = 'c', long = "strict-threads", default_value_t = false)]
    pub strict_threads: bool,

    #[arg(short = 'q', long = "progress", default_value_t = false)]
    pub progress: bool,
}

impl DedupArgs {
    pub fn validate(&self) -> Result<(), Error> {
        if self.percentage < 0.01 || self.percentage > 1.0 {
            anyhow::bail!(
            "Invalid value {} for -p/--percentage! Choose a value between (inclusive) 0.01 and 1.0",
            self.percentage
        )
        }
        Ok(())
    }
}

impl std::fmt::Display for DedupArgs {
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
            self.paired,
            "Percentage".purple(),
            self.percentage,
            "Max edit distance".purple(),
            self.max_edit,
        )?;

        Ok(())
    }
}

const DEDUP_HELP: &str = r#"

RUMINA dedup: cluster and deduplicate or group reads by UMI barcodes

usage:
    rumina dedup -i [*.bam|*.fastq|*.fastq.gz] -g {directional, acyclic, raw} -s <UMI SEPARATOR> [OPTIONS] -o [OUTDIR]

    The input can be one or more BAM or FASTQ files; RUMINA will processall input files sequentially.

arguments:

    [[required]]:
    -i: input files: BAMs must be sorted and indexed.

    -g, --grouping-method: Specifies UMI clustering method. Choose from:
        - directional: as in UMI-tools: predict mutant UMIs based on hamming distance and frequency 
        - acyclic: same as directional, but networks are limited to a depth of one
        - raw: treat UMIs as is; do not error correct

    -s, --separator: Last character in read QNAME immediately before UMI barcode 

    [[clustering]]:
    -l, --length: stratify reads additionally by sequence length (including soft-clipped bases). 
    -u, --rev: search for reverse complements of UMIs when clustering
    -v, --only-group: do not deduplicate clusters; instead annotate reads with cluster ID in BX tag
    -d, --min-depth: minimum number of reads in a cluster for it to be output [3]
    -f, --singletons: remove minimum depth limit for clusters. Identical to --min_depth 1

    [[grouping - advanced]]
    -p, --percentage: The fraction of a parent UMI's read count an offshoot's count must be [3]
    -m, --max-edit: The maximum edit distance delta between two reads for direct linkage 
    
    [[performance, memory]]
    -t, --threads: number of threads to parallelize coordinate processing. Defaults to # sys threads.
    -c, --strict-threads: restrict IO reading to the number specified in --threads. 

    -x, --split-window: Process an input reference alignment by x coordinates at at time.
    Not using this option will process the entire alignment at once. 
    This is recommended when dealing with extreme depth and limited memory.

    -s, --ensure-sorted: only relevant when --split_window is used.
    Ensure output is sorted by buffering all outbound reads until all windows have completed.

    [[paired-end]]
    -l, --paired: Use only R1 for deduplication, and pair output R1 with R2, similar to UMI-tools.

    -m, --merge-pairs: Use a ref fasta to merge overlapping forward/reverse reads with the same UMI.
    Only supports single-reference bams currently. See -b/--min-overlap-bp. 

    -b, --min-overlap-bp: Minimum number of overlapping, matching bases for merging two reads.
    Reads not meeting this criterion will both be discarded. Only relevant with -m/--merge-pairs. 

    [[misc]]
    -o, --outdir: directory (relative to parent of input) in which to store output files.
    -q, --progress: show progress bar, prints to standard error.

"#;
