use crate::args::Args;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::progbars::*;
use crate::read_store::bottomhash::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::Record;
use crate::GroupReport;
use crate::GroupingMethod;
use indicatif::MultiProgress;
use log::info;
use parking_lot::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Debug)]
pub struct Processor {
    pub read_counter: i64,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub only_group: bool,
    pub singletons: bool,
    pub r1_only: bool,
}

impl Processor {
    pub fn new(
        grouping_method: &GroupingMethod,
        group_by_length: bool,
        seed: u64,
        only_group: bool,
        singletons: bool,
        r1_only: bool,
    ) -> Self {
        Processor {
            read_counter: 0,
            min_max: Arc::new(Mutex::new(GroupReport::new())),
            grouping_method: grouping_method.clone(),
            group_by_length,
            seed,
            only_group,
            singletons,
            r1_only,
        }
    }

    pub fn init_from_args(args: &Args, seed: u64) -> Self {
        Self::new(
            &args.grouping_method,
            args.length,
            seed,
            args.only_group,
            args.singletons,
            args.r1_only,
        )
    }

    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads<T: Record + Send + std::fmt::Debug>(
        &mut self,
        bottomhash: &mut BottomHashMap<T>,
        multiprog: &MultiProgress,
        separator: &String,
    ) -> Vec<T> {
        let grouping_method = Arc::new(&self.grouping_method);

        let mut coord_bar = make_coordbar(bottomhash.read_dict.len() as u64);
        coord_bar = multiprog.add(coord_bar);

        coord_bar.set_prefix("REFERENCE COORDINATE");

        let outreads = Arc::new(Mutex::new(Vec::new()));

        bottomhash.read_dict.par_drain(..).for_each(|position| {
            for umi in position.1 {
                let mut umis_reads = umi.1;

                // sort UMIs (stably) by read count in descending order.
                umis_reads
                    .par_sort_by(|_umi1, count1, _umi2, count2| count2.count.cmp(&count1.count));

                let umis = umis_reads
                    .keys()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();

                let umi_len = umis.get(0).unwrap().len();

                let grouper = Grouper::new(&umis, 1, umi_len);
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
                let groupies = grouper.cluster(counts.clone(), Arc::clone(&grouping_method));
                let tagged_reads = group_handler.tag_records(counts, groupies, umis_reads);

                // update grouping report
                let mut out = outreads.lock();
                out.extend(tagged_reads.1);
                drop(out);

                if let Some(group_report) = tagged_reads.0 {
                    let mut min_max = self.min_max.lock();
                    min_max.update(group_report, num_umis);
                    drop(min_max)
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
    // organize reads in bottomhash based on position
    pub fn pull_read<T: Record>(
        &mut self,
        read: T,
        pos: i64,
        key: ReadKey,
        bottomhash: &mut BottomHashMap<T>,
        separator: &String,
    ) {
        bottomhash.update_dict(pos, key.get_key(), read.get_umi(separator), read);
        self.read_counter += 1;
    }
}
