pub mod bam_io;
pub mod bam_reader;
pub mod fastq_io;
pub mod file_io;

pub use crate::io::bam_reader::WindowedBamReader;
pub use crate::io::{
    bam_io::BamIO,
    fastq_io::FastqIO,
    file_io::{gather_files, process_all, FileIO},
};
