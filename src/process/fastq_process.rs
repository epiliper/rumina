use crate::args::Args;
use crate::io::FastqIO;
use crate::io::FileIO;
use crate::process::file_process::FileProcess;
use crate::processor::Processor;
use crate::progbars::ProgressTracker;
use crate::read_store::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::{FastqRecord, Record};
use crate::utils::gen_outfile_name;
use anyhow::{Context, Error};
use colored::Colorize;
use indexmap::IndexMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::Arc;

pub struct FastQFileProcess {
    io: FastqIO,
    chunk_processor: Processor,
    outfile: String,
    separator: String,
    group_reads: bool,
    progress: bool,
}

impl FileProcess for FastQFileProcess {
    fn init_from_args(args: &Args, file_path: &String, file_name: &String) -> Result<Self, Error> {
        let outfile = gen_outfile_name(Some(&args.outdir), ".fastq.gz", "RUMINA", file_name)?;
        let io = FastqIO::init_from_args(args, file_path, &outfile)?;

        let mut hasher = DefaultHasher::new();
        file_name.hash(&mut hasher);
        let seed = hasher.finish();

        let chunk_processor = Processor::init_from_args(args, seed);
        let separator = args.separator.clone();
        let group_reads = args.only_group;
        let progress = args.progress;

        Ok(Self {
            io,
            chunk_processor,
            outfile,
            separator,
            group_reads,
            progress,
        })
    }

    fn process(mut self) -> Result<(), Error> {
        let (mut pos, mut key): (i64, ReadKey);
        let mut outreads: Vec<FastqRecord> = Vec::with_capacity(1_000_000);

        let mut bottomhash: BottomHashMap<FastqRecord> = BottomHashMap {
            read_dict: IndexMap::with_capacity(500),
            read_count: 0,
        };

        let mut reader = self.io.reader.take().context("Reader unitialized")?;
        let mut pt = ProgressTracker::initialize_main(1, self.progress);

        let mut refresh_count = 0;
        pt.initialize_windows(1);
        pt.intake_reads_msg();
        for record in reader.records() {
            let r = record?;
            (pos, key) = r.get_pos_key(self.chunk_processor.group_by_length);
            self.chunk_processor.pull_read(
                r,
                pos,
                key,
                &mut bottomhash,
                &self.separator,
                self.group_reads,
            )?;

            refresh_count += 1;

            if refresh_count == 1000 {
                bottomhash.shrink_to_fit();
                refresh_count = 0;
            }
        }

        pt.update_window_reads(bottomhash.read_count);
        // println! {"Processing {} reads...", bottomhash.read_count};

        assert_eq!(bottomhash.read_dict.keys().len(), 1);
        assert_eq!(
            bottomhash
                .read_dict
                .first()
                .context("Empty read dict")?
                .1
                .len(),
            1
        );

        outreads.extend(self.chunk_processor.group_reads(
            &mut bottomhash,
            &mut pt.coord_bar,
        ));

        self.io.write_reads(&mut outreads);

        let num_reads_in = self.chunk_processor.read_counter;
        let min_maxes = self.chunk_processor.min_max.clone();

        drop(self.chunk_processor);

        // do final report
        let mut group_report = Arc::into_inner(min_maxes)
            .context("Unable to deref group reports")?
            .into_inner();

        group_report.num_reads_input_file = num_reads_in;

        pt.finish();
        // report on min and max number of reads per group
        // this creates minmax.txt
        if !group_report.is_blank() {
            println!("{}", "DONE".green());

            group_report.write_to_report_file(&self.outfile);
            println!("{}\n", group_report);
        }

        Ok(())
    }
}
