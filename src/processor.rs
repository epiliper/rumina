use crate::args::Args;
use crate::deduplicator::GroupHandler;
use crate::grouper::Grouper;
use crate::read_store::bottomhash::BottomHashMap;
use crate::readkey::ReadKey;
use crate::record::Record;
use crate::GroupReport;
use crate::GroupingMethod;
use anyhow::Error;
use indicatif::ProgressBar;
use log::info;
use parking_lot::Mutex;
use rayon::prelude::*;
use smol_str::SmolStr;
use std::collections::HashMap;
use std::sync::Arc;

pub type UmiHistogram<'a> = HashMap<&'a str, (i32, bool)>;

#[derive(Debug)]
pub struct Processor {
    pub read_counter: i64,
    pub min_max: Arc<Mutex<GroupReport>>,
    pub grouping_method: GroupingMethod,
    pub group_by_length: bool,
    pub seed: u64,
    pub only_group: bool,
    pub min_depth: usize,
    pub paired: bool,
    percentage: f32,
    max_edit: u32,
    cluster_rev: bool,
}

impl Processor {
    pub fn new(
        grouping_method: &GroupingMethod,
        group_by_length: bool,
        seed: u64,
        only_group: bool,
        min_depth: usize,
        paired: bool,
        percentage: f32,
        max_edit: u32,
        cluster_rev: bool,
    ) -> Self {
        assert!(percentage > 0.0 && percentage <= 1.0);
        Processor {
            read_counter: 0,
            min_max: Arc::new(Mutex::new(GroupReport::new())),
            grouping_method: grouping_method.clone(),
            group_by_length,
            seed,
            only_group,
            min_depth,
            paired,
            percentage,
            max_edit,
            cluster_rev,
        }
    }

    pub fn init_from_args(args: &Args, seed: u64) -> Self {
        let min_depth;

        if args.singletons {
            min_depth = 1;
        } else {
            min_depth = args.min_cluster_depth;
        }

        Self::new(
            &args.grouping_method,
            args.length,
            seed,
            args.only_group,
            min_depth,
            args.paired,
            args.percentage,
            args.max_edit,
            args.cluster_rev,
        )
    }

    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn group_reads<T: Record + Send + std::fmt::Debug>(
        &mut self,
        bottomhash: &mut BottomHashMap<T>,
        coord_bar: &mut ProgressBar,
    ) -> Vec<T> {
        let grouping_method = Arc::new(&self.grouping_method);

        coord_bar.set_length(bottomhash.read_dict.len() as u64);

        let outreads = Arc::new(Mutex::new(Vec::new()));

        bottomhash
            .read_dict
            .par_drain(..)
            .for_each(|(position, mut key_map)| {
                for (key, mut umi_read_map) in key_map.drain(..) {
                    // sort UMIs (stably) by read count in descending order.
                    umi_read_map.par_sort_by(|_umi1, (count1, _map1), _umi2, (count2, _map2)| {
                        count2.cmp(&count1)
                    });

                    let umis = umi_read_map
                        .keys()
                        .map(|x| x.clone())
                        .collect::<Vec<SmolStr>>();

                    let umi_len = umis.get(0).unwrap().len();

                    let grouper = Grouper::new(
                        &umis,
                        self.max_edit,
                        self.percentage,
                        umi_len,
                        self.cluster_rev,
                    );
                    let mut counts: UmiHistogram = HashMap::with_capacity(umi_read_map.len());

                    // get number of reads for each raw UMI
                    let mut num_umis = 0;

                    for umi in &umis {
                        counts
                            .entry(umi.as_str())
                            .or_insert((umi_read_map[umi].0, true));
                        num_umis += 1;
                    }

                    let mut group_handler = GroupHandler {
                        // make seed for tag unique per position and key
                        seed: self.seed + position as u64 + key,
                        group_only: self.only_group,
                        min_depth: self.min_depth,
                    };

                    // perform UMI clustering per the method specified
                    let groupies = grouper.cluster(counts.clone(), Arc::clone(&grouping_method));
                    let (group_report, tagged_reads) = group_handler
                        .tag_records(groupies, &mut umi_read_map, counts)
                        .unwrap();

                    // update grouping report
                    let mut out = outreads.lock();
                    out.extend(tagged_reads);
                    drop(out);

                    if let Some(group_report) = group_report {
                        let mut min_max = self.min_max.lock();
                        min_max.update(group_report, num_umis);
                        drop(min_max)
                    }
                }
                coord_bar.inc(1);
            });
        // coord_bar.finish_and_clear();

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
        retain_all: bool,
    ) -> Result<(), Error> {
        bottomhash.update_dict(
            pos,
            key.get_key(),
            read.get_umi(separator)?,
            read,
            retain_all,
        );
        self.read_counter += 1;
        Ok(())
    }
}
