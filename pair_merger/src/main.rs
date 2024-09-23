use clap::Parser;
use merge::handle_dupes;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{record::Aux, Format, Header, Read, Reader, Record, Writer};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::string::String;
use std::sync::{
    mpsc::{channel, Receiver, Sender},
    Arc, Mutex,
};
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

    let (tx, rx): (Sender<Record>, Receiver<Record>) = channel();

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("ERROR: Invalid number of threads specified");

    // get list of duplicate UMIs
    let mut duplicate_umis = HashSet::new();

    let dupes = File::open(dupe_list).expect("Duplicate barcode file not found!");
    let dupes = BufReader::new(dupes);

    for line in dupes.lines() {
        let umi = line.unwrap();
        duplicate_umis.insert(umi.trim().to_string());
    }

    // write every 1000 reads to the output bam
    let writer_handle = thread::spawn(move || {
        let mut bam_writer = Writer::from_path(outbam, &header, Format::Bam).unwrap();
        bam_writer.set_threads(args.threads).unwrap();

        let mut buffer = Vec::with_capacity(1000000);
        let mut counter: u32 = 0;

        loop {
            match rx.recv() {
                Ok(read) => {
                    buffer.push(read);
                    counter += 1;

                    if counter == 1000 {
                        buffer
                            .iter()
                            .for_each(|read| bam_writer.write(&read).expect("Error writing read"));
                        buffer.clear();
                        counter = 0;
                    }
                }
                Err(e) => {
                    eprintln!("ERR: {e}");
                    if !buffer.is_empty() {
                        println!("Not empty");
                        buffer
                            .iter()
                            .for_each(|read| bam_writer.write(&read).expect("Error writing read"));
                    }
                    break;
                }
            }
        }
    });

    let holding: Arc<Mutex<HashMap<String, Vec<Record>>>> = Arc::new(Mutex::new(HashMap::new()));

    // check all reads to see if have flagged duplicate UMI
    // if not, write to output
    bam_reader.records().par_bridge().for_each(|r| {
        let read = r.expect("Error reading read in!");

        let umi = if let Ok(Aux::String(bx_i)) = read.aux(UMI_TAG) {
            bx_i
        } else {
            "NULL"
        };

        if duplicate_umis.contains(umi) {
            holding
                .lock()
                .unwrap()
                .entry(umi.to_string())
                .or_insert(Vec::new())
                .push(read);
        } else {
            tx.send(read).expect("Error sending read!");
        }
    });

    let holding = Arc::into_inner(holding)
        .expect("Cannot dereference dupe reads")
        .into_inner()
        .expect("Mutex poisoned!");

    let mut merged_reads = handle_dupes(holding);

    merged_reads
        .drain(0..)
        .for_each(|read| tx.send(read).unwrap());

    drop(tx);
    writer_handle.join().expect("Writer thread panicked!");
}
