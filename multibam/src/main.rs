use crossbeam::channel::{bounded, Receiver, Sender};
use rayon::iter::ParallelDrainRange;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Format, Header, IndexedReader, Read, Record, Writer};
use std::fs;
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
    let args = Args::parse();

    let input = Path::new(&args.input);

    let output_dir = Path::parent(&input).unwrap().join(&args.out);
    if !output_dir.exists() {
        let _ = fs::create_dir(&output_dir);
    }

    let output_dir = fs::canonicalize(&output_dir).unwrap();

    let mut bam_reader = IndexedReader::from_path(&input).unwrap();
    bam_reader.set_threads(args.threads).unwrap();

    let ref_count = bam_reader.header().target_count();
    let window_size = args.split_window as i64;

    // for every reference, split it into windows by split_window
    for tid in 0..ref_count {
        let max_pos = bam_reader.header().target_len(tid).unwrap() as i64;

        let mut ranges: Vec<[i64; 2]> = Vec::new();

        let mut j;
        let mut i = 0;

        while i <= max_pos - window_size {
            j = i;
            i = j + window_size;

            ranges.push([j, i]);
        }

        // add the remainder
        // fetch end value is exclusive, so adding 1
        ranges.push([i, max_pos + 1]);

        // write to a subfile for each window
        ranges.par_drain(..).enumerate().for_each(|(idx, range)| {
            let mut bam_reader = IndexedReader::from_path(&input).unwrap();
            bam_reader.set_threads(args.threads).unwrap();
            let header = Header::from_template(bam_reader.header());

            let mut reads_to_write_f = Vec::with_capacity(1_000_000);
            println!("{:?}", range);

            let start = range[0];
            let end = range[1];

            let _ = bam_reader.fetch((tid, start, end + 300));

            for r in bam_reader.records() {
                let read = r.unwrap();

                if read.is_reverse() {
                    if read.reference_end() <= end && read.reference_end() >= start {
                        reads_to_write_f.push(read);
                    }
                } else if read.reference_start() < end && read.reference_start() >= start {
                    reads_to_write_f.push(read);
                }
            }

            let infile = Path::file_name(input).unwrap().to_str().unwrap();
            let prefix = infile.split(".bam").next().unwrap();

            let outfile_f = format!("{}{}{}{}{}{}", prefix, ".", "f", tid, idx, ".bam");

            let mut bam_writer_f =
                Writer::from_path(&output_dir.join(&outfile_f), &header, Format::Bam).unwrap();

            bam_writer_f.set_threads(args.threads).unwrap();

            reads_to_write_f.iter().for_each(|read| {
                bam_writer_f.write(read).unwrap();
            });
        });
    }
}
