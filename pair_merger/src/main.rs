use clap::Parser;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{record::Aux, Format, Header, Read, Reader, Record, Writer};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::string::String;
use std::sync::mpsc;
use std::sync::{
    mpsc::{Receiver, Sender},
    Arc, Mutex,
};

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

    let (tx, rx): (Sender<Record>, Receiver<Record>) = mpsc::channel();

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("ERROR: Invalid number of threads specified");

    // get list of duplicate UMIs
    let mut duplicate_umis = HashSet::new();

    let dupes = File::open(dupe_list).expect("Duplicate barcode file not found!");

    for line in BufReader::new(dupes).lines() {
        let umi = line.unwrap();
        duplicate_umis.insert(umi.trim().to_string());
    }

    thread::spawn(move || )

    let mut bam_reader = Reader::from_path(inbam).unwrap();
    bam_reader.set_threads(args.threads).unwrap();

    let mut bam_writer = Writer::from_path(
        outbam,
        &Header::from_template(bam_reader.header()),
        Format::Bam,
    )
    .unwrap();
    bam_writer.set_threads(args.threads).unwrap();

    // iterate over all reads in input
    // if they contain dupe tag, then hold them for merging.
    // otherwise, write them to the new file.
    //
    // and yes, use ALL THE THREADS
    //
    let holding: Arc<Mutex<Vec<Record>>> = Arc::new(Mutex::new(Vec::new()));

    bam_reader.records().par_bridge().for_each(|r| {
        let read = r.unwrap();

        let umi = if let Ok(Aux::String(bx_i)) = read.aux(UMI_TAG) {
            bx_i
        } else {
            "NULL"
        };

        if duplicate_umis.contains(umi) {
            holding.lock().unwrap().push(read);
        } else {
            tx.send(read).unwrap();
        }
    });

    todo!()
}
