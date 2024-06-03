/* This module contains scripts for deduplication and error-correction
* of reads within a UMI group.
*
* correct_errors() is the driver function;
* 1. input all the reads within the UMI group
* 2. group the reads by sequence
* 3. select the group with the highest phred score across reads
* 4. output the read with the highest phred score from selected group
*/

extern crate bam;
use crate::bottomhash::ReadsAndCount;
use crate::IndexMap;
use bam::Record;
use std::collections::HashMap;

pub fn correct_errors(clusters: &mut Vec<ReadsAndCount>) -> Record {
    let mut sequences: IndexMap<Box<[u8]>, (Vec<Record>, i32)> = IndexMap::new();
    let _counts: IndexMap<&Box<[u8]>, i32> = IndexMap::new();

    // group the reads by sequence
    for cluster in clusters {
        cluster.reads.drain(0..).for_each(|x| {
            sequences
                .entry(Box::from(*&x.sequence().raw()))
                .or_insert((Vec::new(), 0));
            sequences[x.sequence().raw()].1 += 1;
            sequences[x.sequence().raw()].0.push(x);
        });
    }

    // the sequence with the most reads is at index 0 (can be tied)
    sequences.sort_by(|_a, b, _c, d| d.1.cmp(&b.1));

    let mut phreddies: Vec<Vec<Record>> = Vec::new();
    let mut first = true;
    let mut max = 0;

    // get all sequence groups with the max observed read count (necessary if more than one group
    // has highest num of reads)
    for cluster in sequences.drain(0..) {
        if first {
            max = cluster.1 .1;
            phreddies.push(cluster.1 .0);
            first = false;
        } else if cluster.1 .1 == max {
            phreddies.push(cluster.1 .0);
        } else {
            break;
        }
    }
    // get the group with the best overall phred score
    // and return the best phred score read from it
    return get_best_phred(phreddies);
}

// get the number of reads across all UMIs within a group
// this is useful for setting a threshold for reads observed per UMI group
pub fn get_counts(top_umi: &Vec<&String>, counts: &HashMap<&String, i32>) -> i64 {
    let mut read_count = 0;
    for umi in top_umi {
        read_count += counts.get(umi).unwrap();
    }
    return read_count as i64;
}

// get the read with the best phred score from the read group
pub fn get_best_read(mut reads: Vec<Record>) -> Record {
    let mut quals: IndexMap<i32, Record> = IndexMap::new();
    for read in reads.drain(0..) {
        // need to multiply avg by 100 and convert to int to use hashing
        let avg_phred = ((read
            .qualities()
            .raw()
            .iter()
            .map(|x| *x as f32)
            .sum::<f32>()
            / read.qualities().raw().len() as f32)
            * 100.0) as i32;
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
