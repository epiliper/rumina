use clap::Parser;
use std::path::PathBuf;

pub const PHRED33: core::ops::RangeInclusive<u8> = 33..=77;
pub const SOLEXA: core::ops::RangeInclusive<u8> = 59..=106;
pub const PHRED64: core::ops::RangeInclusive<u8> = 64..=106;

#[derive(Clone, clap::ValueEnum, Debug)]
#[allow(clippy::style)]
pub enum FastqQualEncoding {
    PHRED33,
    SOLEXA,
    PHRED64,
}

impl std::fmt::Display for FastqQualEncoding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FastqQualEncoding::PHRED64 => write!(f, "phred64"),
            FastqQualEncoding::PHRED33 => write!(f, "phred33"),
            FastqQualEncoding::SOLEXA => write!(f, "solexa"),
        }
    }
}

#[derive(Parser, Debug)]
#[clap(author, version, override_help = EXTRACT_HELP)]
pub struct ExtractArgs {
    #[arg(short = 'i')]
    pub in1: PathBuf,

    #[arg(short = 'I')]
    pub in2: Option<PathBuf>,

    #[arg(short = 'o')]
    pub out1: PathBuf,

    #[arg(short = 'O')]
    pub out2: Option<PathBuf>,

    #[arg(short = 's', default_value_t = '_')]
    pub umi_separator: char,

    #[arg(short = 'p')]
    pub pattern1: Option<String>,

    #[arg(short = 'P')]
    pub pattern2: Option<String>,

    #[arg(long = "retain-seq")]
    pub retain_seq: bool,

    #[arg(long = "mask-qual", default_value_t = 0)]
    pub mask_qual: u8,

    #[arg(long = "quality-filter")]
    pub min_qual: Option<u8>,

    #[arg(short = 'e', long = "quality-encoding", default_value_t = FastqQualEncoding::PHRED33)]
    pub qual_encoding: FastqQualEncoding,

    #[arg(short = 'b', long = "batch-size", default_value_t = 1e4 as usize)]
    pub batch_size: usize,

    #[arg(short = 't', long = "threads", default_value_t = num_cpus::get())]
    pub threads: usize,
}

const EXTRACT_HELP: &str = r#"
RUMINA Extract - add UMIs from FASTQ read sequence into read headers

Usage: 
    Single-end:
        rumina extract -i <FASTQ> -p <PATTERN> -o <OUTPUT> [OPTIONS]

    Paired-end:
        rumina extract -i <FASTQ> -p <PATTERN> -o <OUTPUT1>
            -I <FASTQ2> -P <PATTERN2> -O <OUTPUT2>

    All FASTQ files must end with either .fastq, or .fastq.gz.

    If both R1 and R2 input FASTQs are provided, you must also specify a pattern for each (-p and -P).

Options:
    misc:
        --version: show version

    extract:
        -i: first FASTQ input. Must end in .fastq or .fastq.gz
        -I: second FASTQ input. Must end in .fastq or .fastq.gz
        -o: first FASTQ output. Must end in .fastq or .fastq.gz
        -O: second FASTQ output. Must end in .fastq. or .fastq.gz

        -p: extraction pattern for file given with -i
        -P: extraction pattern for file given with -I.

        -s: character to use to delimit barcodes from each other and read header ['_']

        --retain-seq: don't remove barcode bases from read sequences during extraction.
        Barcode sequence will be in both the read header and sequence.

        --mask-qual: replace any UMI barcode bases below this quality with 'N'

        --quality-filter: don't output reads with UMIs with base(s) below this quality.
        In paired-end data, if one mate fails this filter, the other will be removed.

        -e/--qual-encoding: quality encoding to use for filtering/masking. 
        Choose from "phred33", "phred64", or "solexa".

        -b/--batch-size: number of reads to buffer before writing [10000]

        -t/--threads: number of threads to use for compression/parallel extraction. [Number of system threads]
"#;
