use crate::args::Args;
use crate::io::BamIO;
use crate::io::FileIO;
use crate::pair_merger::PairMerger;
use crate::process::file_process::FileProcess;
use crate::processor::Processor;
use crate::progbars::ProgressTracker;
use crate::read_store::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::{BamRecord, Record};
use crate::utils::{gen_outfile_name, index_bam};
use anyhow::Error;
use colored::Colorize;
use indexmap::IndexMap;
use log::info;
use std::fs::remove_file;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::Arc;

pub struct BamFileProcess {
    io: BamIO,
    chunk_processor: Processor,
    outfile: String,
    pair_merger: Option<PairMerger>,
    separator: String,
    group_reads: bool,
}

impl FileProcess for BamFileProcess {
    fn init_from_args(args: &Args, file_path: &String, file_name: &String) -> Result<Self, Error> {
        let outfile = gen_outfile_name(Some(&args.outdir), ".bam", "RUMINA", file_name)?;
        let bam_io = BamIO::init_from_args(args, file_path, &outfile);

        let mut hasher = DefaultHasher::new();
        file_name.hash(&mut hasher);
        let seed = hasher.finish();

        let chunk_processor = Processor::init_from_args(args, seed);
        let mut pair_merger: Option<PairMerger> = None;
        let group_reads = args.only_group;

        if let Some(ref ref_fasta) = args.merge_pairs {
            pair_merger = Some(PairMerger {
                ref_fasta: ref_fasta.to_string(),
                min_overlap_bp: args.min_overlap_bp,
                threads: args.threads,
                infile: outfile.to_string(),
                outfile: gen_outfile_name(None, ".bam", "MERGED", &outfile)?,
                split_window: args.split_window,
            })
        }

        let separator = args.separator.clone();

        Ok(Self {
            io: bam_io,
            chunk_processor,
            outfile,
            pair_merger,
            separator,
            group_reads,
        })
    }

    fn process(mut self) -> Result<(), Error> {
        let (mut pos, mut key): (i64, ReadKey);
        let mut outreads: Vec<BamRecord> = Vec::with_capacity(1_000_000);

        let mut pt =
            ProgressTracker::initialize_main(self.io.windowed_reader.meta_header.target_count());

        while self.io.windowed_reader.next_reference()? {
            pt.initialize_windows(self.io.windowed_reader.windows.len());
            while self.io.windowed_reader.next_window() {
                let mut bottomhash = BottomHashMap {
                    read_dict: IndexMap::with_capacity(500),
                    read_count: 0,
                };

                let mut window_records = 0;
                pt.intake_reads_msg();

                for record in self.io.windowed_reader.window_records() {
                    if !record.is_first_in_template() && self.chunk_processor.r1_only {
                        continue;
                    }

                    (pos, key) = record.get_pos_key(self.chunk_processor.group_by_length);
                    self.chunk_processor.pull_read(
                        record,
                        pos,
                        key,
                        &mut bottomhash,
                        &self.separator,
                        self.group_reads,
                    )?;
                    window_records += 1;
                }

                info!("{} reads pulled from window", window_records);
                pt.update_window_reads(window_records);

                outreads.extend(Processor::group_reads(
                    &mut self.chunk_processor,
                    &mut bottomhash,
                    &mut pt.coord_bar,
                    &self.separator,
                ));

                self.io.write_reads(&mut outreads);

                info!(
                    "Processed {} total reads...",
                    self.chunk_processor.read_counter
                );
                pt.next_window()
            }
            pt.next_ref();
        }

        pt.finish();

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

        drop(self.io.writer); // dropping to avoid vague samtools warning
        let idx = index_bam(&self.outfile, self.io.num_threads).expect("Failed to index bam");

        if let Some(mut pair_merger) = self.pair_merger {
            info!("{:?}", pair_merger);

            let merge_report = pair_merger.merge_windows();
            remove_file(self.outfile).ok();
            remove_file(idx).ok();
            index_bam(&pair_merger.outfile, self.io.num_threads).unwrap();
            print!("{merge_report}");
        }

        Ok(())
    }
}
