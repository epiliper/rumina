use bam::bam_writer::BamWriterBuilder;
use bam::{Record, RecordWriter};
use indexmap::IndexMap;
use rayon::prelude::*;
use std::fs;
use std::ops::Range;
use std::path::Path;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(term_width = 10)]
struct Args {
    #[arg(long = "in", index = 1)]
    input: String,

    #[arg(long = "out", index = 2)]
    out: String,

    #[arg(long = "split_window", index = 3)]
    split_window: usize,

    #[arg(long = "threads", index = 4)]
    threads: usize,
}

fn main() {
    // let args: Vec<String> = env::args().collect();
    //
    let args = Args::parse();

    let input = Path::new(&args.input);
    let add_threads = (args.threads - 1) as u16;

    let output_dir = Path::parent(&input).unwrap().join(&args.out);
    if !output_dir.exists() {
        fs::create_dir(&output_dir);
    }

    let output_dir = fs::canonicalize(&output_dir).unwrap();

    let mut input_bam = bam::BamReader::from_path(&input, add_threads).unwrap();
    let input_header = bam::BamReader::from_path(&input, 0)
        .unwrap()
        .header()
        .clone();
    let max_pos = input_header.reference_len(0).unwrap() as i32;
    let mut pos_range = (0..max_pos).collect::<Vec<i32>>();

    // let chunk_size = &args[3].parse::<usize>().unwrap();
    let chunk_size = args.split_window;

    let mut ranges_of_reads: IndexMap<Range<i32>, Vec<Record>> = IndexMap::new();
    let mut ranges: Vec<Range<i32>> = Vec::new();

    for idx_array in pos_range.chunks_mut(chunk_size) {
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

    ranges_of_reads
        .par_values_mut()
        .enumerate()
        .for_each(|(idx, read_bundle)| {
            let infile = Path::file_name(input).unwrap().to_str().unwrap();
            let prefix = infile.split(".bam").next().unwrap();
            let outfile = format!("{}{}{}{}", prefix, ".", idx, ".bam");

            print! {"\rWriting {:?} reads to {:?}", read_bundle.len(), outfile};

            let mut output_writer = BamWriterBuilder::from_path(
                &mut BamWriterBuilder::new().additional_threads(add_threads),
                &output_dir.join(outfile),
                input_header.clone(),
            )
            .unwrap();

            read_bundle
                .drain(0..)
                .for_each(|x| output_writer.write(&x).unwrap());
        });
}
