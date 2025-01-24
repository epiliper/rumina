use crate::group_report::StaticUMI;
use crate::readkey::ReadKey;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::index::{build, Type};
use rust_htslib::bam::{Header, IndexedReader, Read, Record, Writer};
use std::fs::DirEntry;
use std::path::Path;

pub fn get_umi<'b>(record: &'b Record, separator: &String) -> &'b str {
    unsafe {
        std::str::from_utf8_unchecked(record.qname())
            .rsplit_once(separator)
            .expect("ERROR: failed to get UMI from read QNAME. Check --separator. Exiting.")
            .1
    }
}

pub fn get_file_ext(file: &DirEntry) -> String {
    file.file_name()
        .to_str()
        .and_then(|s| s.rsplit_once('.'))
        .map(|(_, ext)| ext.to_string())
        .unwrap_or_default()
}

pub fn get_umi_static<'c>(raw_umi: &'c str) -> StaticUMI {
    let mut umi = StaticUMI::new();
    umi.extend(raw_umi.as_bytes().into_iter().map(|b| *b));
    umi
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

pub fn get_windows(
    window_size: Option<i64>,
    bam: &IndexedReader,
    reference_id: u32,
) -> Vec<[i64; 2]> {
    if let Some(window_size) = window_size {
        if window_size == 0 {
            return vec![[0, i64::MAX]];
        }

        let max_pos = bam
            .header()
            .target_len(reference_id)
            .expect("Invalid reference") as i64;

        let mut ranges: Vec<[i64; 2]> = Vec::new();

        let mut j;
        let mut i = 0;

        while i <= max_pos - window_size {
            j = i;
            i = j + window_size;

            ranges.push([j, i]);
        }

        ranges.push([i, max_pos + 1]);

        return ranges;
    } else {
        return vec![[0, i64::MAX]];
    }
}

pub fn get_read_pos_key(group_by_length: bool, read: &Record) -> (i64, ReadKey) {
    let mut pos;
    let key: ReadKey;

    if read.is_reverse() {
        // set end pos as start to group with forward-reads covering same region
        pos = read.reference_end();
        pos += read.cigar().trailing_softclips(); // pad with right-side soft clip
        key = ReadKey {
            length: read.seq_len() * group_by_length as usize,
            reverse: true,
            chr: read.tid() as usize,
        };
        (pos, key)
    } else {
        pos = read.reference_start();
        pos -= read.cigar().leading_softclips(); // pad with left-side soft clip

        key = ReadKey {
            length: read.seq_len() * group_by_length as usize,
            reverse: false,
            chr: read.tid() as usize,
        };
        (pos, key)
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

pub fn index_bam(bam_name: &String, num_threads: usize) -> Result<(), rust_htslib::errors::Error> {
    build(
        Path::new(bam_name),
        None,
        Type::Bai,
        num_threads.try_into().unwrap(),
    )
}
