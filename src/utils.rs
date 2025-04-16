use crate::group_report::StaticUMI;
use rust_htslib::bam::{index, Header, IndexedReader, Read, Writer};
use std::path::Path;

pub fn get_umi_static<'c>(raw_umi: &'c str) -> StaticUMI {
    let mut umi = StaticUMI::new();
    umi.extend(raw_umi.as_bytes().into_iter().map(|b| *b));
    umi
}

pub fn get_file_ext(path: &Path) -> Option<&str> {
    path.extension().and_then(|ext| ext.to_str())
}

pub fn gen_outfile_name(outdir: Option<&String>, suffix: &str, fname: &str) -> String {
    let outf = fname.rsplit_once(".bam").unwrap().0.to_string() + &format!("_{suffix}.bam");
    if let Some(outdir) = outdir {
        Path::new(outdir)
            .join(Path::new(&outf))
            .to_str()
            .unwrap()
            .to_string()
    } else {
        Path::new(&outf).to_str().unwrap().to_string()
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

        return ranges;
    } else {
        return vec![Window {
            start: 0,
            end: i64::MAX,
        }];
    }
}

pub fn make_bam_writer(file_name: &String, header: Header, num_threads: usize) -> Writer {
    let mut bam_writer =
        Writer::from_path(file_name, &header, rust_htslib::bam::Format::Bam).unwrap();
    bam_writer.set_threads(num_threads).unwrap();
    bam_writer
}

pub fn make_bam_reader(input_file: &String, num_threads: usize) -> (Header, IndexedReader) {
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
        Ok(_) => return Ok(idx_name),
        Err(e) => return Err(e),
    }
}
