use crate::args::Args;
use crate::bam_io::{bam_io::BamIO, fastq_io::FastqIO};
use crate::main_dedup::process_windows;
use crate::pair_merger::PairMerger;
use crate::utils::index_bam;
use crate::utils::{gen_outfile_name, get_file_ext};
use crate::window_processor::Processor;
use colored::Colorize;
use log::{error, info};
use std::collections::HashMap;
use std::fs::{read_dir, remove_file};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use std::sync::Arc;

pub fn gather_files(input_file: &str) -> HashMap<String, String> {
    let inpath = Path::new(input_file);

    if inpath.is_dir() {
        read_dir(inpath)
            .into_iter()
            .flatten()
            .filter_map(|entry| {
                let entry = match entry {
                    Ok(e) => e,
                    Err(e) => {
                        error!("Unable to read file: {}. Won't be used in processing.", e);
                        return None;
                    }
                };
                let path = entry.path();
                if !path.is_dir() && get_file_ext(&path) == Some("bam") {
                    Some((
                        path.to_string_lossy().into_owned(),
                        entry.file_name().to_string_lossy().into_owned(),
                    ))
                } else {
                    None
                }
            })
            .collect()
    } else {
        std::iter::once((
            inpath.to_string_lossy().into_owned(),
            inpath
                .file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_default(),
        ))
        .collect()
    }
}

pub trait FileProcess {
    fn init_from_args(args: &Args, bam_file_path: &String, bam_file_name: &String) -> Self;
    fn process(self);
}

pub struct BamFileProcess {
    io: BamIO,
    chunk_processor: Processor,
    outfile: String,
    pair_merger: Option<PairMerger>,
    separator: String,
}

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

    fn process(self) {
        // unimplemented!();
        info!("{:?}", self.chunk_processor);
    }
}

impl FileProcess for BamFileProcess {
    fn init_from_args(args: &Args, file_path: &String, file_name: &String) -> Self {
        let outfile = gen_outfile_name(Some(&args.outdir), "RUMINA", file_name);
        let bam_io = BamIO::init_from_args(args, file_path, &outfile);

        let mut hasher = DefaultHasher::new();
        file_name.hash(&mut hasher);
        let seed = hasher.finish();

        let chunk_processor = Processor::init_from_args(args, seed);
        let mut pair_merger: Option<PairMerger> = None;

        if let Some(ref ref_fasta) = args.merge_pairs {
            pair_merger = Some(PairMerger {
                ref_fasta: ref_fasta.to_string(),
                min_overlap_bp: args.min_overlap_bp,
                threads: args.threads,
                infile: outfile.to_string(),
                outfile: gen_outfile_name(None, "MERGED", &outfile),
                split_window: args.split_window,
            })
        }

        let separator = args.separator.clone();

        Self {
            io: bam_io,
            chunk_processor,
            outfile,
            pair_merger,
            separator,
        }
    }

    fn process(mut self) {
        info!("{:?}", self.chunk_processor);

        process_windows(
            &mut self.chunk_processor,
            self.io.windowed_reader,
            self.io.mate_reader,
            &self.separator,
            self.io.writer,
        );

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

        let idx = index_bam(&self.outfile, self.io.num_threads).expect("Failed to index bam");

        if let Some(mut pair_merger) = self.pair_merger {
            info!("{:?}", pair_merger);

            let merge_report = pair_merger.merge_windows();
            remove_file(self.outfile).ok();
            remove_file(idx).ok();
            index_bam(&pair_merger.outfile, self.io.num_threads).unwrap();
            print!("{merge_report}");
        }
    }
}
