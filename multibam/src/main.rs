use indexmap::IndexMap;
use rust_htslib::bam::{ext::BamRecordExtensions, Format, Header, Read, Reader, Record, Writer};
use std::env;
use std::fs;
use std::ops::Range;
use std::path::Path;

fn main() {
    let args: Vec<String> = env::args().collect();
    let input = Path::new(&args[1]);
    let output_dir = Path::parent(input).unwrap().join(&args[2]);
    if !output_dir.exists() {
        fs::create_dir(&output_dir);
    }

    let output_dir = fs::canonicalize(&output_dir).unwrap();

    let mut input_bam = Reader::from_path(&args[1]).unwrap();
    let header = Header::from_template(&input_bam.header());
    let input_header = header.to_hashmap();

    let max_pos = input_header.get("SQ").unwrap()[0].get("LN").unwrap();
    let mut pos_range = (0..max_pos.parse::<i64>().unwrap()).collect::<Vec<i64>>();

    let chunk_size = &args[3].parse::<usize>().unwrap();

    let mut ranges_of_reads: IndexMap<Range<i64>, Vec<Record>> = IndexMap::new();
    let mut ranges: Vec<Range<i64>> = Vec::new();

    for idx_array in pos_range.chunks_mut(*chunk_size) {
        let read_range = Range {
            start: *idx_array.iter().min().unwrap(),
            end: *idx_array.iter().max().unwrap(),
        };

        print! {"\r{:?}", read_range};

        ranges_of_reads
            .entry(read_range.clone())
            .or_insert(Vec::new());
        ranges.push(read_range);
    }

    let _ = input_bam.set_threads(4);

    for read in input_bam.records() {
        let mut pos = i64::MAX;
        let r = read.unwrap();
        if r.is_reverse() {
            pos = r.reference_end();
        } else {
            pos = r.reference_start();
        }

        for range in &ranges {
            match range.contains(&pos) {
                true => {
                    ranges_of_reads[range].push(r);
                    break;
                }
                false => {}
            };
        }
    }

    for (i, set) in ranges_of_reads.values_mut().enumerate() {
        let infile = Path::file_name(input).unwrap().to_str().unwrap();
        let prefix = infile.split(".bam").next().unwrap();
        let outfile = format!("{}{}{}{}", prefix, ".", i, ".bam");

        print! {"\rWriting {:?} reads to {:?}", set.len(), outfile};

        let mut writer =
            Writer::from_path(&output_dir.join(outfile), &header, Format::Bam).unwrap();
        let _ = writer.set_threads(4);

        set.drain(0..).for_each(|x| writer.write(&x).unwrap());
    }
}
