use crate::utils::get_drain_end;
use indexmap::IndexMap;
use rayon::slice::ParallelSliceMut;
use rust_htslib::bam::{Record, Writer};

pub struct ReadChannel {
    reads: IndexMap<i64, Vec<Record>>,
    min_pos: i64,
    writer: Writer,
}

impl ReadChannel {
    pub fn new(writer: Writer) -> ReadChannel {
        ReadChannel {
            reads: IndexMap::new(),
            min_pos: 0,
            writer,
        }
    }

    pub fn write_sorted_to_file(&mut self, write_all: bool) {
        // to avoid breaking sorting, write only forward reads that start BEFORE the reverse the
        // reads
        // this will inevitably leave some forward reads behind, which should be written
        // afterwards.
        //
        //
        let mut outreads: Vec<Record> = Vec::new();

        if write_all {
            for read in self.reads.drain(..).flat_map(|(_pos, vec)| vec) {
                outreads.push(read);
            }
        } else {
            self.update();

            self.reads.par_sort_keys();

            // let end_idx = get_drain_end(&self.reads, self.min_pos + 1);
            let end_idx = self.reads.len() / 2;
            println!("{}", end_idx);

            for (_pos, f_read) in self.reads.drain(..end_idx) {
                outreads.extend(f_read);
            }
        }

        outreads.par_sort_by(|ra, rb| ra.pos().cmp(&rb.pos()));
        println!("Writing {} reads", outreads.len());
        for read in outreads {
            self.writer.write(&read).unwrap();
        }
        self.update();
    }

    pub fn update(&mut self) {
        self.min_pos = *self.reads.keys().min().unwrap_or(&0);
    }

    pub fn add_read(&mut self, read: Record) {
        self.reads
            .entry(read.pos())
            .or_insert_with(|| Vec::new())
            .push(read);
    }

    pub fn intake(&mut self, reads: Vec<Record>, write_all: bool) {
        reads.into_iter().for_each(|read| self.add_read(read));

        self.write_sorted_to_file(write_all);
    }
}
