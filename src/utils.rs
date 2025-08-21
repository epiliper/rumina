use anyhow::{Context, Error};
use rust_htslib::bam::{index, Header, IndexedReader, Read, Writer};
use std::path::Path;

pub struct RecordFile {
    pub fname: String,
    pub fpath: String,
}

pub enum FileType {
    BamFile(RecordFile),
    FastqFile(RecordFile),
}

pub fn identify_file_type(path: &Path) -> Option<FileType> {
    let fname = path.file_name()?.to_str()?.to_string();
    let fpath = path.to_str()?.to_string();

    if fname.ends_with(".fastq.gz") {
        return Some(FileType::FastqFile(RecordFile { fname, fpath }));
    }

    if fname.ends_with(".bam") {
        return Some(FileType::BamFile(RecordFile { fname, fpath }));
    }

    None
}

pub fn gen_outfile_name(
    outdir: Option<&String>,
    split: &str,
    suffix: &str,
    fname: &str,
) -> Result<String, Error> {
    let outf = fname
        .rsplit_once(split)
        .context("Unable to split file")?
        .0
        .to_string()
        + &format!("_{suffix}{split}");
    if let Some(outdir) = outdir {
        Ok(Path::new(outdir)
            .join(Path::new(&outf))
            .to_str()
            .context("unable to construct output file path")?
            .to_string())
    } else {
        Ok(Path::new(&outf)
            .to_str()
            .context("Unable to construct output file path")?
            .to_string())
    }
}

#[derive(Clone, Debug)]
pub struct Window {
    pub start: i64,
    pub end: i64,
}

pub fn get_windows(window_size: Option<i64>, max_pos: i64) -> Vec<Window> {
    if let Some(window_size) = window_size {
        if window_size == 0 {
            return vec![Window {
                start: 0,
                end: i64::MAX,
            }];
        }

        let mut ranges: Vec<Window> = Vec::new();

        let mut j;
        let mut i = 0;

        while i <= max_pos - window_size {
            j = i;
            i = j + window_size;

            ranges.push(Window { start: j, end: i });
        }

        ranges.push(Window {
            start: i,
            end: max_pos + 1,
        });

        ranges
    } else {
        vec![Window {
            start: 0,
            end: i64::MAX,
        }]
    }
}

pub fn make_bam_writer(file_name: &str, header: Header, num_threads: usize) -> Writer {
    let mut bam_writer =
        Writer::from_path(file_name, &header, rust_htslib::bam::Format::Bam).unwrap();
    bam_writer.set_threads(num_threads).unwrap();
    bam_writer
}

pub fn make_bam_reader(input_file: &str, num_threads: usize) -> (Header, IndexedReader) {
    let mut bam_reader = IndexedReader::from_path(input_file).unwrap();
    bam_reader.set_threads(num_threads).unwrap();
    let header = Header::from_template(bam_reader.header());

    (header, bam_reader)
}

pub fn index_bam(
    bam_name: &String,
    num_threads: usize,
) -> Result<String, rust_htslib::errors::Error> {
    // this function will return an error if the input bam is not sorted.
    let idx_name = format!("{bam_name}.bai");
    let res = index::build(
        Path::new(bam_name),
        Some(Path::new(&idx_name)),
        index::Type::Bai,
        num_threads.try_into().unwrap(),
    );

    match res {
        Ok(_) => Ok(idx_name),
        Err(e) => Err(e),
    }
}
