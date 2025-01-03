use crate::bottomhash::BottomHashMap;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::progbars::*;
use crate::readkey::ReadKey;
use crate::report::BarcodeSender;
use crate::report::BarcodeTracker;
use crate::utils::get_umi;
use crate::GroupReport;
use crate::GroupingMethod;
use indicatif::MultiProgress;
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::{Record, Writer};
use std::collections::HashMap;
use std::sync::Arc;

pub struct ChunkProcessor {
    pub read_counter: i64,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub split_window: Option<i64>,
    pub only_group: bool,
    pub singletons: bool,
    pub track_barcodes: Option<String>,
    pub r1_only: bool,
    pub barcode_tracker: Arc<Mutex<BarcodeTracker>>,
}

impl ChunkProcessor {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads(
        &mut self,
        bottomhash: &mut BottomHashMap,
        bc_sender: &Option<BarcodeSender>,
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
                    track_barcodes: self.track_barcodes.is_some(),
                };

                // perform UMI clustering per the method specified
                let groupies = grouper.cluster(counts, Arc::clone(&grouping_method));

                let tagged_reads = group_handler.tag_records(
                    groupies,
                    umis_reads,
                    Arc::clone(&self.barcode_tracker),
                );

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

            if let Some(s) = bc_sender {
                for (bc, count) in self.barcode_tracker.lock().barcode_counter.drain(..) {
                    s.send((bc, count))
                        .expect("failed to send barcode to writing channel!");
                }
            }
            coord_bar.inc(1);
        });

        coord_bar.finish_and_clear();
        Arc::try_unwrap(outreads)
            .expect("Unable to dereference tagged reads!")
            .into_inner()
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

    pub fn write_reads(&mut self, mut outreads: Vec<Record>, bam_writer: &mut Writer) {
        outreads
            .drain(..)
            .for_each(|read| bam_writer.write(&read).unwrap());
    }
}
