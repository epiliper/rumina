use crate::args::Args;
use crate::cli::print_file_info;
use crate::process::{BamFileProcess, FastQFileProcess, FileProcess};
use crate::record::Record;
use crate::utils::get_file_ext;
use anyhow::{Context, Error, Result};
use std::collections::HashMap;
use std::fs::read_dir;
use std::path::Path;

pub trait FileIO<T: Record> {
    fn write_reads(&mut self, outreads: &mut Vec<T>);
}

pub fn gather_files(input_file: &str) -> Result<HashMap<String, String>, anyhow::Error> {
    let inpath = Path::new(input_file);

    if inpath.is_dir() {
        let entries = read_dir(inpath)
            .with_context(|| format!("Failed to read directory: {}", input_file))?;

        Ok(entries
            .into_iter()
            .filter_map(|entry| {
                let entry = entry.ok()?;
                let path = entry.path();

                if !path.is_dir() {
                    let ext = get_file_ext(&path);

                    if ext == Some("bam") || ext == Some("gz") {
                        return Some((
                            path.to_string_lossy().into_owned(),
                            entry.file_name().to_string_lossy().into_owned(),
                        ));
                    } else {
                        println! {"Skipping file {:?}; unrecognized extension", entry.file_name()};
                        return None;
                    }
                } else {
                    return None;
                }
            })
            .collect())
    } else {
        if !inpath.exists() {
            return Err(anyhow::anyhow!("File does not exist: {}", input_file));
        }

        Ok(std::iter::once((
            inpath.to_string_lossy().into_owned(),
            inpath
                .file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .ok_or_else(|| anyhow::anyhow!("Failed to get file name for: {}", input_file))?,
        ))
        .collect())
    }
}

// Assuming get_file_ext is defined elsewhere, ensure it doesn't panic
// pub fn process_all(args: &Args, file_map: HashMap<String, String>) -> Result<(), Error> {
pub fn process_all(args: &Args, file_map: HashMap<String, String>) -> Vec<Result<(), Error>> {
    let num_files = file_map.len();
    let mut results = vec![];

    for (i, (file_path, file_name)) in file_map.into_iter().enumerate() {
        print_file_info(&file_name, i + 1, num_files);
        let p = match get_file_ext(Path::new(&file_name)) {
            Some("bam") => BamFileProcess::init_from_args(args, &file_path, &file_name)
                .and_then(|process| process.process()),
            Some("fastq.gz") => FastQFileProcess::init_from_args(args, &file_path, &file_name)
                .and_then(|process| process.process()),
            _ => Ok(()),
        };

        results.push(p.with_context(|| format!("Failed to process file {}", file_name)));
    }

    results
}
