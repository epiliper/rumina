use arrayvec::ArrayVec;
use colored::Colorize;
use indexmap::IndexMap;
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::ops::AddAssign;
use std::path::Path;

const MAX_UMI_LENGTH: usize = 30;

// This module defines the report generated during deduplication.
#[derive(Debug)]
pub struct BarcodeTracker {
    barcode_counter: IndexMap<StaticUMI, u16>,
}

// a UMI barcode in an ArrayVec of u8, efficient for
// reporting barcodes without operating on them.
// #[derive(Debug)]
pub type StaticUMI = ArrayVec<u8, MAX_UMI_LENGTH>;

impl BarcodeTracker {
    pub fn new() -> Self {
        BarcodeTracker {
            barcode_counter: IndexMap::with_capacity(1000),
        }
    }

    // count UMI
    pub fn count(&mut self, umi: StaticUMI) {
        self.barcode_counter.entry(umi).or_insert(0).add_assign(1);
    }

    // update UMI counts after last round of deduplication
    pub fn update(&mut self, other_barcode_tracker: &mut BarcodeTracker) {
        for (umi, count) in self.barcode_counter.drain(..) {
            other_barcode_tracker
                .barcode_counter
                .entry(umi)
                .and_modify(|existing_count| *existing_count += count)
                .or_insert(count);
        }
    }
}

// This report contains details like UMIs in/out, reads in/out, and other details, and is updated
// after the deduplication of each batch.
#[derive(Debug)]
pub struct GroupReport {
    pub min_reads: i64,
    pub max_reads: i64,
    pub min_reads_group: [u8; 8],
    pub max_reads_group: [u8; 8],
    pub num_passing_groups: i64,
    pub num_groups: i64,
    pub num_umis: i64,
    pub num_reads_input_file: i64,
    pub num_reads_output_file: i64,
    pub barcode_tracker: BarcodeTracker,
}

impl GroupReport {
    pub fn new() -> Self {
        GroupReport {
            min_reads: i64::MAX,
            min_reads_group: *b"NONENONE",
            max_reads: 0,
            max_reads_group: *b"NONENONE",
            num_passing_groups: 0,
            num_groups: 0,
            num_umis: 0,
            num_reads_input_file: 0,
            num_reads_output_file: 0,
            barcode_tracker: BarcodeTracker::new(),
        }
    }

    // detect if report is empty (occurs if no groups pass singleton filter)
    pub fn is_blank(&self) -> bool {
        self.max_reads == 0
    }

    // after a batch has been processed, check to see if fields need to be udpated
    pub fn update(&mut self, mut other_report: GroupReport, num_umis: i32) {
        if other_report.max_reads > self.max_reads {
            self.max_reads = other_report.max_reads;
            self.max_reads_group = other_report.max_reads_group;
        }

        if other_report.min_reads < self.min_reads {
            self.min_reads = other_report.min_reads;
            self.min_reads_group = other_report.min_reads_group;
        }

        // count the number of UMI groups used in consensus
        self.num_passing_groups += other_report.num_passing_groups;
        self.num_groups += other_report.num_groups;
        self.num_umis += num_umis as i64;

        // record the number of reads to be written
        self.num_reads_output_file += other_report.num_reads_output_file;
        // self.barcode_tracker.update(other_report.barcode_tracker);
        other_report
            .barcode_tracker
            .update(&mut self.barcode_tracker);
    }

    // once deduplication of the file is complete, only list UMIs that were observed more than
    // once.
    pub fn write_to_report_file(&mut self, output_file: &String, barcode_file: bool) {
        let report_file = Path::new(&output_file).parent().unwrap().join("minmax.txt");
        //
        if !report_file.exists() {
            let _ = File::create(&report_file);
        }

        let mut report_f = OpenOptions::new()
            .append(true)
            .open(report_file)
            .expect("unable to open minmax file");

        let _ = report_f.write(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                self.num_reads_input_file,
                self.num_reads_output_file,
                self.num_umis,
                self.num_groups,
                self.num_passing_groups,
                String::from_utf8(self.min_reads_group.to_vec()).unwrap(),
                self.min_reads,
                String::from_utf8(self.max_reads_group.to_vec()).unwrap(),
                self.max_reads,
            )
            .as_bytes(),
        );

        if barcode_file {
            let barcode_file = Path::new(&output_file)
                .parent()
                .unwrap()
                .join("barcodes.tsv");

            if !barcode_file.exists() {
                let _ = File::create(&barcode_file);
            }

            let mut barcode_f = OpenOptions::new()
                .append(true)
                .open(barcode_file)
                .expect("unable to open minmax file");

            for (umi, count) in self.barcode_tracker.barcode_counter.iter() {
                writeln!(barcode_f, "{}\t{}", String::from_utf8_lossy(&umi), count).unwrap();
            }
        }
    }
}

// printed after file completion
impl fmt::Display for GroupReport {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}",
            "Minimum reads per group".cyan(),
            self.min_reads,
            "Maximum reads per group".cyan(),
            self.max_reads,
            "Total UMI groups".cyan(),
            self.num_groups,
            "Groups passing singleton filtering".cyan(),
            self.num_passing_groups,
            "Total UMIs considered".cyan(),
            self.num_umis,
            "Input reads (mapped)".cyan(),
            self.num_reads_input_file,
            "Output reads".cyan(),
            self.num_reads_output_file
        )
    }
}
