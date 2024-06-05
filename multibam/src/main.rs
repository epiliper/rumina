use indexmap::IndexMap;
use bam::bam_writer::BamWriterBuilder;
use bam::bgzip::Writer;
use bam::{Record, RecordWriter};
use std::env;
use std::fs;
use std::ops::Range;
use std::path::Path;
use rayon::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();
    let input = Path::new(&args[1]);
    let output_dir = Path::parent(input).unwrap().join(&args[2]);
    if !output_dir.exists() {
        fs::create_dir(&output_dir);
    }

    let output_dir = fs::canonicalize(&output_dir).unwrap();

    let mut input_bam = bam::BamReader::from_path(&args[1], 4).unwrap();
    let input_header = bam::BamReader::from_path(&args[1], 0).unwrap().header().clone();
    let max_pos = input_header.reference_len(0).unwrap() as i32;
    // let input_header = header.to_hashmap();

    // let max_pos = input_header.get("SQ").unwrap()[0].get("LN").unwrap();
    // let mut pos_range = (0..max_pos.parse::<i64>().unwrap()).collect::<Vec<i64>>();
    let mut pos_range = (0..max_pos).collect::<Vec<i32>>();

    let chunk_size = &args[3].parse::<usize>().unwrap();

    let mut ranges_of_reads: IndexMap<Range<i32>, Vec<Record>> = IndexMap::new();
    let mut ranges: Vec<Range<i32>> = Vec::new();

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

    for read in input_bam {
        let mut pos = i32::MAX;
        let r = read.unwrap();
        if r.flag().is_reverse_strand() {
            pos = r.calculate_end();
        } else {
            pos = r.start();
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

    ranges_of_reads.par_values_mut().enumerate().for_each(
        
        |(idx, read_bundle)| {
            let infile = Path::file_name(input).unwrap().to_str().unwrap();
            let prefix = infile.split(".bam").next().unwrap();
            let outfile = format!("{}{}{}{}", prefix, ".", idx, ".bam");

            print! {"\rWriting {:?} reads to {:?}", read_bundle.len(), outfile};

            let mut output_writer = BamWriterBuilder::from_path(&mut BamWriterBuilder::new().additional_threads(8),
            &output_dir.join(outfile), 
            input_header.clone()).unwrap();

            read_bundle.drain(0..).for_each(|x| output_writer.write(&x).unwrap());
        }

    );

}
