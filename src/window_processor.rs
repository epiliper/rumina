use crate::bottomhash::BottomHashMap;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::progbars::*;
use crate::readkey::ReadKey;
use crate::utils::{get_umi, Window};
use crate::GroupReport;
use crate::GroupingMethod;
use indexmap::IndexSet;
use indicatif::MultiProgress;
use log::info;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read, Record, Writer};
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Debug)]
pub struct ChunkProcessor {
    pub read_counter: i64,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub split_window: Option<i64>,
    pub only_group: bool,
    pub singletons: bool,
    pub r1_only: bool,
}

impl ChunkProcessor {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(
        &mut self,
        bottomhash: &mut BottomHashMap,
        multiprog: &MultiProgress,
        separator: &String,
    ) -> Vec<Record> {
        let grouping_method = Arc::new(&self.grouping_method);

        let mut coord_bar = make_coordbar(bottomhash.read_dict.len() as u64);
        coord_bar = multiprog.add(coord_bar);

        coord_bar.set_prefix("REFERENCE COORDINATE");

        let outreads = Arc::new(Mutex::new(Vec::new()));

        bottomhash.read_dict.par_drain(..).for_each(|position| {
            for umi in position.1 {
                let mut umis_reads = umi.1;

                // sort UMIs by read count
                // note that this is an unstable sort, so we need to identify read-tied groups
                umis_reads
                    .par_sort_by(|_umi1, count1, _umi2, count2| count2.count.cmp(&count1.count));

                let umis = umis_reads
                    .keys()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();

                let grouper = Grouper { umis: &umis };
                let mut counts: HashMap<&str, i32> = HashMap::with_capacity(umis_reads.len());

                // get number of reads for each raw UMI
                let mut num_umis = 0;

                for umi in &umis {
                    counts.entry(umi.as_str()).or_insert(umis_reads[umi].count);
                    num_umis += 1;
                }

                let mut group_handler = GroupHandler {
                    seed: self.seed + position.0 as u64, // make seed unique per position
                    group_only: self.only_group,
                    singletons: self.singletons,
                    separator: &separator,
                };

                // perform UMI clustering per the method specified
                let groupies = grouper.cluster(counts, Arc::clone(&grouping_method));

                let tagged_reads = group_handler.tag_records(groupies, umis_reads);

                drop(grouper);

                // update grouping report
                if let Some(tagged_reads) = tagged_reads {
                    let mut out = outreads.lock();
                    out.extend(tagged_reads.1);
                    drop(out);

                    if let Some(group_report) = tagged_reads.0 {
                        let mut min_max = self.min_max.lock();
                        min_max.update(group_report, num_umis);
                        drop(min_max)
                    }
                }
            }

            coord_bar.inc(1);
        });

        coord_bar.finish_and_clear();
        info!("Outputting final reads for writing...");
        info!("\n{:?}", self.min_max);

        Arc::try_unwrap(outreads)
            .expect("Unable to dereference tagged reads!")
            .into_inner()
    }

    pub fn retrieve_r2s(
        tid: u32,
        chunk_start: i64,
        chunk_end: i64,
        reader: &mut IndexedReader,
        ids: IndexSet<&[u8]>,
    ) -> Vec<Record> {
        reader.fetch((tid, chunk_start, chunk_end + 15)).unwrap();
        let mut mates: Vec<Record> = Vec::with_capacity(ids.len());

        for read in reader.records().flatten() {
            if ids.contains(read.qname()) && read.is_last_in_template() && read.pos() >= chunk_start
            {
                mates.push(read);
            }
        }

        mates
    }

    // organize reads in bottomhash based on position
    pub fn pull_read(
        &mut self,
        read: Record,
        pos: i64,
        key: ReadKey,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        bottomhash.update_dict(
            pos,
            key.get_key(),
            get_umi(&read, separator).to_string(),
            read,
        );
        self.read_counter += 1;
    }

    pub fn write_reads(
        &mut self,
        outreads: &mut Vec<Record>,
        bam_writer: &mut Writer,
        bam_reader: &mut Option<IndexedReader>,
        tid: u32,
        window: &Window,
    ) {
        let mut count = 0;
        let mut mates: Option<Vec<Record>> = None;

        if !outreads.is_empty() {
            if let Some(ref mut bam_reader) = bam_reader {
                let chunk_start = window.start;
                let chunk_end = window.end;

                let mut ids_to_pair = IndexSet::with_capacity(outreads.len());

                outreads.iter().for_each(|read| {
                    ids_to_pair.insert(read.qname());
                });

                mates = Some(Self::retrieve_r2s(
                    tid,
                    chunk_start,
                    chunk_end,
                    bam_reader,
                    ids_to_pair,
                ));
            }

            if let Some(mates) = mates {
                outreads.extend(mates);
            }

            outreads.par_sort_by(|ra, rb| ra.pos().cmp(&rb.pos()));
            outreads.drain(..).for_each(|read| {
                bam_writer.write(&read).unwrap();
                count += 1;
            });
        }

        info!("Written {count} reads!")
    }
}
