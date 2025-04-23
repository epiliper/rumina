use crate::args::Args;
use crate::cli::print_file_info;
use crate::process::{BamFileProcess, FastQFileProcess, FileProcess};
use crate::record::Record;
use crate::utils::{identify_file_type, FileType};
use anyhow::{Context, Error, Result};
use std::fs::read_dir;
use std::path::Path;

pub trait FileIO<T: Record> {
    fn write_reads(&mut self, outreads: &mut Vec<T>);
}

pub fn gather_files(input_file: &str) -> Result<Vec<FileType>, anyhow::Error> {
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
                    if let Some(ftype) = identify_file_type(&path) {
                        Some(ftype)
                    } else {
                        println! {"Skipping file {:?}; unrecognized extension", entry.file_name()};
                        None
                    }
                } else {
                    println! {"Skipping folder {:?}", entry.file_name()};
                    return None;
                }
            })
            .collect())
    } else {
        if !inpath.exists() {
            return Err(anyhow::anyhow!("File does not exist: {}", input_file));
        }

        Ok(std::iter::once(
            identify_file_type(&inpath)
                .ok_or_else(|| anyhow::anyhow!("Unrecognized file extension: {}", input_file)),
        )
        .flatten()
        .collect())
    }
}

pub fn process_all(args: &Args, file_map: Vec<FileType>) -> Vec<Result<(), Error>> {
    let num_files = file_map.len();
    let mut results = vec![];

    for (i, file) in file_map.into_iter().enumerate() {
        let p = match file {
            FileType::BamFile(f) => {
                print_file_info(&f.fname, i + 1, num_files);
                BamFileProcess::init_from_args(args, &f.fpath, &f.fname)
                    .and_then(|process| process.process())
                    .with_context(|| format!("Failed to process file {}", &f.fname))
            }

            FileType::FastqFile(f) => {
                print_file_info(&f.fname, i + 1, num_files);
                FastQFileProcess::init_from_args(args, &f.fpath, &f.fname)
                    .and_then(|process| process.process())
                    .with_context(|| format!("Failed to process file {}", &f.fname))
            }
        };

        results.push(p);
    }

    results
}
