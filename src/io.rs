pub mod bam_io;
pub mod bam_reader;

pub mod fastq_dedup_io;
pub mod fastq_extract_io;

pub mod fastq_output;
pub mod fastq_writer;
pub mod fastqio;

pub mod file_io;

pub use bam_reader::WindowedBamReader;
pub use file_io::FileIO;
