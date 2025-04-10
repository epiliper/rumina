use crate::args::Args;
use crate::bam_io::bam_io::BamIO;
use crate::main_dedup::process_chunks;
use crate::pair_merger::PairMerger;
use crate::utils::index_bam;
use crate::utils::{gen_outfile_name, get_file_ext};
use crate::window_processor::ChunkProcessor;
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

pub struct BamFileProcess {
    bam_io: BamIO,
    chunk_processor: ChunkProcessor,
    // pair_merger: Option<PairMerger>,
}

impl BamFileProcess {
    pub fn init_from_args(args: &Args, outfile: &String, seed: u64) -> Self {
        let bam_io = BamIO::init_from_args(args, outfile);
        let chunk_processor = ChunkProcessor::init_from_args(args, seed);

        Self {
            bam_io,
            chunk_processor,
        }
    }
}

pub fn process(input_file: (String, String), args: &Args) {
    let output_file = gen_outfile_name(Some(&args.outdir), "RUMINA", &input_file.1);

    let mut hasher = DefaultHasher::new();
    input_file.hash(&mut hasher);
    let seed = hasher.finish();

    let mut file_process = BamFileProcess::init_from_args(args, &output_file, seed);

    info!("{:?}", file_process.chunk_processor);

    process_chunks(
        &mut file_process.chunk_processor,
        file_process.bam_io.windowed_reader,
        file_process.bam_io.mate_reader,
        &args.separator,
        file_process.bam_io.writer,
    );

    let num_reads_in = file_process.chunk_processor.read_counter;
    let min_maxes = file_process.chunk_processor.min_max.clone();

    drop(file_process.chunk_processor);

    // do final report
    let mut group_report = Arc::try_unwrap(min_maxes).unwrap().into_inner();
    group_report.num_reads_input_file = num_reads_in;

    // report on min and max number of reads per group
    // this creates minmax.txt
    if !group_report.is_blank() {
        println!("{}", "DONE".green());

        group_report.write_to_report_file(&output_file);
        println!("{}", group_report);
    }

    let idx = index_bam(&output_file, args.threads).expect("Failed to index bam");

    if let Some(ref ref_fasta) = args.merge_pairs {
        let mut merger = PairMerger {
            ref_fasta: ref_fasta.to_string(),
            min_overlap_bp: args.min_overlap_bp,
            threads: args.threads,
            infile: output_file.to_string(),
            outfile: gen_outfile_name(None, "MERGED", &output_file),
            split_window: args.split_window,
        };

        info!("{:?}", merger);

        let merge_report = merger.merge_windows();
        remove_file(output_file).ok();
        remove_file(idx).ok();
        index_bam(&merger.outfile, args.threads).unwrap();
        print!("{merge_report}");
    }
}
