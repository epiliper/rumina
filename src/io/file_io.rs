use crate::args::Args;
use crate::cli::print_file_info;
use crate::process::{BamFileProcess, FastQFileProcess, FileProcess};
use crate::record::Record;
use crate::utils::get_file_ext;
use log::error;
use std::collections::HashMap;
use std::fs::read_dir;
use std::path::Path;

pub trait FileIO<T: Record> {
    fn write_reads(&mut self, outreads: &mut Vec<T>);
}

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
                if !path.is_dir()
                    && (get_file_ext(&path) == Some("bam")
                        || get_file_ext(&path) == Some("fastq.gz"))
                {
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

pub fn process_all(args: &Args, file_map: HashMap<String, String>) {
    let num_files = file_map.len();
    for (i, (file_path, file_name)) in file_map.into_iter().enumerate() {
        print_file_info(&file_name, i + 1, num_files);
        match get_file_ext(Path::new(&file_name)) {
            Some("bam") => BamFileProcess::init_from_args(args, &file_path, &file_name).process(),
            Some("fastq.gz") => {
                FastQFileProcess::init_from_args(args, &file_path, &file_name).process()
            }
            _ => (),
        };
    }
}
