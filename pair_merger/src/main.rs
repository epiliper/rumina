#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use clap::Parser;
use merge::handle_dupes;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{record::Aux, Format, Header, Read, Reader, Record, Writer};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::string::String;
use std::sync::{Arc, Mutex};

use crossbeam::channel::{bounded, Receiver, Sender};

use std::thread;

mod merge;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short = 'i')]
    inbam: String,
    #[arg(short = 'o')]
    outbam: String,
    #[arg(short = 'l')]
    dupe_list: String,
    #[arg(long = "threads")]
    threads: usize,
}

const UMI_TAG: &[u8; 2] = b"BX";

fn main() {
    let args = Args::parse();
    let inbam = args.inbam;
    let outbam = args.outbam;
    let dupe_list = args.dupe_list;

    let mut bam_reader = Reader::from_path(inbam).unwrap();
    let header = Header::from_template(bam_reader.header());

    bam_reader.set_threads(args.threads).unwrap();

    let (s, r): (Sender<Record>, Receiver<Record>) = bounded(10000);

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("ERROR: Invalid number of threads specified");

    // get list of duplicate UMIs
    let mut duplicate_umis: HashMap<String, i32> = HashMap::new();

    let dupes = File::open(dupe_list).expect("Duplicate barcode file not found!");
    let dupes = BufReader::new(dupes);

    for line in dupes.lines() {
        let umi = line.unwrap();
        let (umi, count) = umi
            .split_once("\t")
            .expect("Error: barcode list improperly formatted!");

        let count = count.parse::<i32>().expect("invalid count");

        duplicate_umis
            .entry(umi.trim().to_string())
            .and_modify(|curr_count| *curr_count += count)
            .or_insert(count);
    }

    duplicate_umis.retain(|_umi, count| *count > 1);

    // write every 1000 reads to the output bam
    let writer_handle = thread::spawn(move || {
        let mut bam_writer = Writer::from_path(outbam, &header, Format::Bam).unwrap();
        bam_writer.set_threads(args.threads).unwrap();

        let mut buffer = Vec::with_capacity(1000000);
        let mut counter: u32 = 0;
        let mut num_writes = 0;

        loop {
            match r.recv() {
                Ok(read) => {
                    buffer.push(read);
                    counter += 1;
                    if counter == 1000 {
                        for read in &buffer {
                            bam_writer.write(read).expect("Error writing read");
                            num_writes += 1;
                        }
                        buffer.clear();
                        counter = 0;
                    }
                }
                Err(_) => {
                    if !buffer.is_empty() {
                        for read in &buffer {
                            bam_writer.write(read).expect("Error writing read");
                            num_writes += 1;
                        }
                    }
                    println!("Pair merger: Written {} reads", num_writes);
                    break;
                }
            }
        }
    });

    let holding: Arc<Mutex<HashMap<String, Vec<Record>>>> = Arc::new(Mutex::new(HashMap::new()));

    // check all reads to see if have flagged duplicate UMI
    // if not, write to output
    bam_reader.records().for_each(|r| {
        let read = r.expect("Error reading read in!");

        let umi = if let Ok(Aux::String(bx_i)) = read.aux(UMI_TAG) {
            bx_i
        } else {
            "NULL"
        };

        if duplicate_umis.contains_key(umi) {
            holding
                .lock()
                .expect("Unable to lock!")
                .entry(umi.to_string())
                .or_insert(Vec::new())
                .push(read);
        } else {
            s.send(read).expect("Error sending read");
        }
    });

    let mut remainder = Arc::into_inner(holding)
        .expect("Cannot dereference dupe reads")
        .into_inner()
        .expect("Mutex poisoned!");

    let mut merged_reads = handle_dupes(&mut remainder);
    println!("Merged {} forward-reverse pairs.", merged_reads.len());

    merged_reads
        .drain(..)
        .for_each(|read| s.send(read).expect("Error sending read!"));

    drop(s);

    writer_handle.join().expect("Writer thread panicked!");
}
