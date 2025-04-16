pub mod bam_process;
pub mod fastq_process;
pub mod file_process;

pub use crate::process::{
    bam_process::BamFileProcess, fastq_process::FastQFileProcess, file_process::FileProcess,
};
