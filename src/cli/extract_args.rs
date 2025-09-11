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
