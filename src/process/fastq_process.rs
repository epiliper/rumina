use crate::args::Args;
use crate::io::FastqIO;
use crate::io::FileIO;
use crate::process::file_process::FileProcess;
use crate::processor::Processor;
use crate::read_store::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::{FastqRecord, Record};
use crate::utils::gen_outfile_name;
use colored::Colorize;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::Arc;

pub struct FastQFileProcess {
    io: FastqIO,
    chunk_processor: Processor,
    outfile: String,
    separator: String,
}

impl FileProcess for FastQFileProcess {
    fn init_from_args(args: &Args, file_path: &String, file_name: &String) -> Self {
        let outfile = gen_outfile_name(Some(&args.outdir), "RUMINA", file_name);
        let io = FastqIO::init_from_args(args, file_path, &outfile);

        let mut hasher = DefaultHasher::new();
        file_name.hash(&mut hasher);
        let seed = hasher.finish();

        let chunk_processor = Processor::init_from_args(args, seed);
        let separator = args.separator.clone();

        Self {
            io,
            chunk_processor,
            outfile,
            separator,
        }
    }

    fn process(mut self) {
        let (mut pos, mut key): (i64, ReadKey);
        let mut outreads: Vec<FastqRecord> = Vec::with_capacity(1_000_000);

        let mut bottomhash: BottomHashMap<FastqRecord> = BottomHashMap {
            read_dict: IndexMap::with_capacity(500),
        };

        let reader = self.io.reader.take().expect("Reader unitialized");

        for record in reader.records().flatten() {
            (pos, key) = record.get_pos_key(self.chunk_processor.group_by_length);
            self.chunk_processor
                .pull_read(record, pos, key, &mut bottomhash, &self.separator);
        }

        outreads.extend(self.chunk_processor.group_reads(
            &mut bottomhash,
            &MultiProgress::new(),
            &self.separator,
        ));

        self.io.write_reads(&mut outreads);

        let num_reads_in = self.chunk_processor.read_counter;
        let min_maxes = self.chunk_processor.min_max.clone();

        drop(self.chunk_processor);

        // do final report
        let mut group_report = Arc::try_unwrap(min_maxes).unwrap().into_inner();
        group_report.num_reads_input_file = num_reads_in;

        // report on min and max number of reads per group
        // this creates minmax.txt
        if !group_report.is_blank() {
            println!("{}", "DONE".green());

            group_report.write_to_report_file(&self.outfile);
            println!("{}\n", group_report);
        }
    }
}
