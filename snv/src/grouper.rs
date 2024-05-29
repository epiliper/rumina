extern crate bam;
use crate::bottomhash::ReadsAndCount;
use crate::IndexMap;
use bam::record::cigar::Operation;
use bam::Record;
use polars::functions::concat_df_horizontal;
use polars::lazy::dsl::when;
use polars::prelude::*;
use rand::rngs::ThreadRng;
use rand::Rng;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use std::collections::HashMap;
use std::collections::HashSet;

const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ\
                            abcdefghijklmnopqrstuvwxyz\
                            0123456789)(*&^%$#@!~";

const UMI_TAG_LEN: usize = 8;

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper {}

pub fn generate_tag(
    rng: &mut ThreadRng,
    used_tags: &mut HashSet<[u8; UMI_TAG_LEN]>,
) -> [u8; UMI_TAG_LEN] {
    let ug_tag: [u8; UMI_TAG_LEN] = (0..UMI_TAG_LEN)
        .map(|_| {
            let idx = rng.gen_range(0..CHARSET.len());
            CHARSET[idx]
        })
        .collect::<Vec<u8>>()
        .try_into()
        .unwrap();
    if !used_tags.contains(&ug_tag) {
        return ug_tag;
    } else {
        generate_tag(rng, used_tags)
    }
}

pub fn report_read(read: &Record, group: &String) -> DataFrame {
    // get read length and mapping quality
    let l = vec![read.cigar().calculate_query_len()];
    let group = vec![group.clone()];
    let q = vec![read.mapq() as u32];

    let mut min_qual = u8::MAX;
    let mut max_qual = 0;
    let mut mean_qual = 0;

    read.qualities().raw().iter().for_each(|x| {
        mean_qual += x;
        min_qual = *x.min(&min_qual);
        max_qual = *x.max(&max_qual);
    });

    let mean_qual = (mean_qual as f32) / (read.qualities().raw().len() as f32);
    let max_qual = max_qual as u32;
    let min_qual = min_qual as u32;

    // initialize other dataframe columns
    let mut num_matches = 0;
    let mut soft_clip = 0;
    let mut num_insertions = 0;
    let mut num_deletions = 0;
    let mut num_mismatches = 0;
    let mut hard_clip = 0;
    let mut mismatch_quals: Vec<u32> = Vec::new();
    let mut mismatch_pos: Vec<u32> = Vec::new();

    // read cigar for M to get number of nt aligned to reference
    read.cigar().iter().for_each(|m| match m.1 {

        Operation::SeqMatch => num_matches += m.0,

        Operation::Soft => soft_clip += m.0,

        Operation::Hard => hard_clip += m.0,

        Operation::Insertion => {
            num_insertions += m.0;
            let idx = (num_matches + num_mismatches - num_deletions + hard_clip + soft_clip + num_insertions - 1) as usize;

            mismatch_quals.push(
                read.qualities().raw()[idx]
                    as u32,
            );
            mismatch_pos.push(idx.try_into().unwrap());
        }

        Operation::Deletion => {
            num_deletions += m.0;
        }

        Operation::SeqMismatch => {
            num_mismatches += m.0;
            let idx = (num_matches + num_mismatches - num_deletions + hard_clip + soft_clip + num_insertions - 1) as usize;
            mismatch_quals.push(
                *read.qualities().raw().get(idx).expect(&format!("{}\n {}     {}      {}      {}      {}      {}\n idx: {}\n length of sequence = {}", read.cigar().to_string(), num_matches, num_mismatches, num_deletions, num_insertions, hard_clip, soft_clip, idx, l.get(0).unwrap())) as u32
            );
            mismatch_pos.push(idx.try_into().unwrap());
        }

        _ => (),
    });

    let mut cols = vec![
        Series::new("length", l),
        Series::new("mapq", q),
        Series::new("group", group),
        Series::new("min_phred", [min_qual]),
        Series::new("mean_phred", [mean_qual]),
        Series::new("max_phred", [max_qual]),
        Series::new("matches", [num_matches]),
        Series::new("soft_clip", [soft_clip]),
        Series::new("mismatches", [num_mismatches]),
        Series::new("deletions", [num_deletions]),
        Series::new("insertions", [num_insertions]),
        Series::new("mismatch_quals", mismatch_quals),
        Series::new("mismatch_positions", mismatch_pos),
    ];

    let dfs = cols
        .drain(0..)
        .map(|x| DataFrame::new(vec![x]).unwrap())
        .collect::<Vec<DataFrame>>();

    // make per read dataframe
    let df = concat_df_horizontal(&dfs).unwrap();
    return df;
}

// get the number of reads across all UMIs within a group
pub fn get_counts(top_umi: &Vec<&String>, counts: &HashMap<&String, i32>) -> i64 {
    let mut read_count = 0;
    for umi in top_umi {
        read_count += counts.get(umi).unwrap();
    }

    return read_count as i64;
}

pub fn get_best_read(mut reads: Vec<Record>) -> Record {

    let mut quals: IndexMap<i32, Record> = IndexMap::new();

    for read in reads.drain(0..) {
        let avg_phred = ((read.qualities().raw().iter().map(|x| *x as f32).sum::<f32>() / read.qualities().raw().len() as f32) * 100.0) as i32;
        quals.insert(avg_phred, read);
    }
     
    let winning_read_index = *quals.iter_mut().map(|x| x.0).max().unwrap();
    return quals.swap_remove(&winning_read_index).unwrap();

}

pub fn get_best_phred(mut clusters: Vec<Vec<Record>>) -> Record {
    // return a single read, since this is per-position and we expect one read per UMI group
    match clusters.len() {
        1 => return clusters.drain(0..).next().unwrap().remove(0),

        _ => {
            let mut mean_phreds: IndexMap<i32, Vec<Record>> = IndexMap::new();

            for cluster in clusters.drain(0..) {
                let mut avgs: Vec<f32> = Vec::new();

                for read in &cluster {
                    avgs.push(
                        read.qualities()
                            .raw()
                            .iter()
                            .map(|x| *x as f32)
                            .sum::<f32>()
                            / read.qualities().raw().len() as f32,
                    );
                }

                let cluster_avg = ((avgs.iter().sum::<f32>() / avgs.len() as f32) * 100.0) as i32;
                mean_phreds.insert(cluster_avg, cluster);
            }

            // get the index of the sequences with the best phred score across all reads
            let x = *mean_phreds.iter_mut().map(|x| x.0).max().unwrap();

            let mut best_phred = mean_phreds.swap_remove(&x).unwrap();

            if best_phred.len() == 1 {
                return best_phred.remove(0);
            } else {

                return get_best_read(best_phred);

            }
            
            return mean_phreds.swap_remove(&x).unwrap().remove(0);
        }
    }
}

pub fn correct_errors(clusters: &mut Vec<ReadsAndCount>) -> Record {
    let mut sequences: IndexMap<Box<[u8]>, (Vec<Record>, i32)> = IndexMap::new();
    let mut counts: IndexMap<&Box<[u8]>, i32> = IndexMap::new();

    // group the reads by sequence
    for mut cluster in clusters {
        cluster.reads.drain(0..).for_each(|x| {
            sequences
                .entry(Box::from(*&x.sequence().raw()))
                .or_insert((Vec::new(), 0));
            sequences[x.sequence().raw()].1 += 1;
            sequences[x.sequence().raw()].0.push(x);
        });
    }

    // the sequence with the most reads is at index 0 (can be tied)
    sequences.sort_by(|a, b, c, d| d.1.cmp(&b.1));

    let mut phreddies: Vec<Vec<Record>> = Vec::new();

    let mut first = true;
    let mut max = 0;
    for cluster in sequences.drain(0..) {

        if first {
            max = cluster.1.1;
            phreddies.push(cluster.1.0);
            first = false;
        } else if cluster.1.1 == max {
            phreddies.push(cluster.1.0);
        } else {

            break;
        }

    }

    // count number of reads per unique sequence
    // for read_group in & mut sequences {
    //     counts
    //         .entry(*&read_group.0)
    //         .or_insert(read_group.1 .1);
    // }

    // let max = counts.values().next().unwrap();

    // // gather all groups that have maximum reported per-group read count
    // let phreddies: Vec<Vec<Record>> = Vec::new();

    // let mut groups_to_sort: HashSet<&Box<[u8]>> = HashSet::new();
    // for cluster in sequences.drain(0..) {
    //     if cluster.1 .1 == *max {
    //         phreddies.push(sequences.shift_remove(&cluster.0).unwrap().0);
    //         // i'm cloning a box so this should work
    //         groups_to_sort.insert(&cluster.0);
    //     } else {
    //         break;
    //     }
    // }

    // // remove read clusters and send for phred comparison
    // let phreddies = groups_to_sort
    //     .iter()
    //     .map(|x| sequences.shift_remove(*x).unwrap().0)
    //     .collect::<Vec<Vec<Record>>>();
    return get_best_phred(phreddies);
}

impl Grouper {
    // remove the reads associated with each UMI from the bundle
    // tag them
    // push them to a list of tagged Records awaiting writing to an output bamfile
    pub fn tag_groups(
        &mut self,
        final_umis: Vec<Vec<&String>>,
        umis_records: &mut IndexMap<String, ReadsAndCount>,
        counts: HashMap<&String, i32>,
    ) -> ((Vec<DataFrame>, Option<IndexMap<i64, i64>>), Vec<Record>) {
        // for each UMI within a group, assign the same tag
        let mut report_list: Vec<DataFrame> = Vec::new();
        // let mut min_max_report = (i64::MAX, 0);

        let mut rng = rand::thread_rng();
        let mut output_list: Vec<Record> = Vec::with_capacity(1_000_000);
        let mut used_tags: HashSet<[u8; UMI_TAG_LEN]> = HashSet::with_capacity(output_list.len());

        let mut first = true;

        let mut min_max_report: IndexMap<i64, i64> = IndexMap::new();

        for top_umi in final_umis {
            let num_reads_in_group = get_counts(&top_umi, &counts);
            if num_reads_in_group >= 3 {

                if first {
                    first = false;
                    min_max_report.insert(0, i64::MAX);
                    min_max_report.insert(1, 0);

                }
                // check if number of reads per group is new minimum or maximum
                if num_reads_in_group < *min_max_report.get(&0).unwrap() {
                    min_max_report.entry(0).and_modify(|x| *x = num_reads_in_group);
                } 

                if num_reads_in_group > *min_max_report.get(&1).unwrap() {
                    min_max_report.entry(1).and_modify(|x| *x = num_reads_in_group);
                }

                let mut cluster_list: Vec<ReadsAndCount> = Vec::new();
                let ug_tag = generate_tag(&mut rng, &mut used_tags);
                for group in &top_umi {
                    cluster_list.push(umis_records.swap_remove(*group).unwrap());
                }

                    // this is the read from the read group with the highest average phred score
                let mut final_record = correct_errors(&mut cluster_list);

                // final_records.drain(0..).for_each(|mut x| {
                //     // gather data
                //     // report_list.push(report_read(&x, group));
                //     x.tags_mut().push_string(b"UG", &ug_tag);
                //     x.tags_mut().push_string(b"BX", group.as_bytes());
                //     output_list.push(x);
                // })

                final_record.tags_mut().push_string(b"UG", &ug_tag);
                final_record.tags_mut().push_string(b"BX", &top_umi.iter().next().unwrap().as_bytes());
                output_list.push(final_record);
                // }
            }
        }

        if min_max_report.is_empty() {
            
            return ((report_list, None), output_list);
        }

        return ((report_list, Some(min_max_report)), output_list);
    }

    // driver function of grouping
    // recieves outout of main_grouper()
    pub fn tag_records(
        &mut self,
        grouping_output: (HashMap<&String, i32>, Option<Vec<Vec<&String>>>),
        mut umis_records: IndexMap<String, ReadsAndCount>,
    ) -> Option<((Vec<DataFrame>, Option<IndexMap<i64, i64>>), Vec<Record>)> {

        match grouping_output.1 {
            Some(groups) => {
                let reads = self.tag_groups(groups, &mut umis_records, grouping_output.0);
                return Some(reads);
            }
            None => return None,
        }
    }
}
