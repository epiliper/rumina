/* This module contains scripts for deduplication and error-correction
* of reads within a UMI group.
*
* correct_errors() is the driver function;
* 1. input all the reads within the UMI group
* 2. group the reads by sequence
* 3. select the group with the highest phred score across reads
* 4. Output one read from the group
*/

use crate::bottomhash::ReadsAndCount;
use crate::IndexMap;
use rust_htslib::bam::Record;
use std::borrow::Cow;
use std::collections::HashMap;

pub fn correct_errors(clusters: &mut Vec<ReadsAndCount>) -> Vec<Record> {
    if clusters.len() == 1 && clusters[0].reads.len() == 1 {
        return vec![clusters.remove(0).reads.remove(0)];
    }

    let mut sequences: IndexMap<Cow<'static, [u8]>, (Vec<Record>, i32)> =
        IndexMap::with_capacity(clusters.len());

    let _counts: IndexMap<&Box<[u8]>, i32> = IndexMap::new();

    // group the reads by sequence
    for cluster in clusters {
        cluster.reads.drain(..).for_each(|x| {
            sequences
                .entry(Cow::Owned(x.seq().encoded.to_vec()))
                .or_insert((Vec::new(), 0));
            sequences[x.seq().encoded].1 += 1;
            sequences[x.seq().encoded].0.push(x);
        });
    }

    // the sequence with the most reads is at index 0 (can be tied)
    sequences.sort_by(|a, b, c, d| d.1.cmp(&b.1).then_with(|| a.cmp(c)));

    let mut most_reads_groups: Vec<Vec<Record>> = Vec::new();
    let mut first = true;
    let mut max = 0;

    // get all sequence groups with the max observed read count (necessary if more than one group
    // has highest num of reads)
    for cluster in sequences.drain(..) {
        if first {
            max = cluster.1 .1;
            most_reads_groups.push(cluster.1 .0);
            first = false;
        } else if cluster.1 .1 == max {
            most_reads_groups.push(cluster.1 .0);
        // once we hit groups with fewer reads, stop searching
        } else {
            break;
        }
    }
    // get the group with the best overall phred score
    // and return the best phred score read from it
    vec![get_best_phred(most_reads_groups)]
}

// used with the --group_only arg to return all reads within a group with a group tag
pub fn push_all_reads(clusters: &mut Vec<ReadsAndCount>) -> Vec<Record> {
    let mut reads_to_write: Vec<Record> = Vec::with_capacity(clusters.len());

    clusters
        .drain(..)
        .for_each(|read_group| reads_to_write.extend(read_group.reads));

    reads_to_write
}

// get the number of reads across all UMIs within a group
// this is useful for setting a threshold for reads observed per UMI group
pub fn get_counts(top_umi: &Vec<&str>, counts: &HashMap<&str, i32>) -> i64 {
    let mut read_count = 0;
    for umi in top_umi {
        read_count += counts.get(umi).unwrap();
    }
    read_count as i64
}

pub fn get_best_phred(mut clusters: Vec<Vec<Record>>) -> Record {
    // return a single read, since this is per-position and we expect one read per UMI group
    match clusters.len() {
        // if just one group, return one read
        1 => {
            let mut cluster = clusters.drain(..).next().unwrap();
            cluster.sort_by(|ra, rb| ra.qname().cmp(rb.qname()));
            return cluster.remove(0);
        }

        // if two groups are tied for read majority, pick the one with the best overall phred score
        _ => {
            let mut mean_phreds: IndexMap<i32, Vec<Record>> =
                IndexMap::with_capacity(clusters.len());

            for cluster in clusters.drain(..) {
                let mut avgs: Vec<f32> = Vec::new();

                for read in &cluster {
                    avgs.push(
                        read.qual().iter().map(|x| *x as f32).sum::<f32>()
                            / read.qual().len() as f32,
                    );
                }

                let cluster_avg = ((avgs.iter().sum::<f32>() / avgs.len() as f32) * 100.0) as i32;
                mean_phreds.insert(cluster_avg, cluster);
            }

            // mean_phreds.par_sort_keys();
            let x = *mean_phreds.iter_mut().map(|x| x.0).max().unwrap();

            let mut best_phred = mean_phreds.swap_remove(&x).unwrap();

            // remove one read from this group (final read to represent UMI group)
            best_phred.sort_by(|ra, rb| rb.qname().cmp(&ra.qname()));
            best_phred.swap_remove(0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bottomhash::ReadsAndCount;
    use rust_htslib::bam::Record;
    use std::collections::HashMap;

    #[test]
    fn test_correct_errors_single_read() {
        let mut clusters: Vec<ReadsAndCount> = Vec::new();
        let mut reads = Vec::new();
        let mut record = Record::new();
        // Simulate a BAM record with some sequence and quality scores
        record.set(b"read1", None, b"ATCG", b"####");
        reads.push(record);
        clusters.push(ReadsAndCount { reads, count: 1 });
        let result = correct_errors(&mut clusters);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_correct_errors_multiple_reads() {
        let mut clusters: Vec<ReadsAndCount> = Vec::new();
        let mut reads1 = Vec::new();
        let mut reads2 = Vec::new();
        let mut record1 = Record::new();
        let mut record2 = Record::new();
        let mut record3 = Record::new();
        // Simulate BAM records with some sequences and quality scores
        record1.set(b"read1", None, b"ATCG", b"####");
        record2.set(b"read2", None, b"ATCG", b"####");
        record3.set(b"read3", None, b"ATGG", b"####");
        reads1.push(record1);
        reads1.push(record2);
        reads2.push(record3);
        clusters.push(ReadsAndCount {
            reads: reads1,
            count: 2,
        });
        clusters.push(ReadsAndCount {
            reads: reads2,
            count: 1,
        });
        let result = correct_errors(&mut clusters);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_get_counts() {
        let top_umi = vec!["UMI1", "UMI2"];
        let counts = HashMap::from([("UMI1", 10), ("UMI2", 5)]);
        let result = get_counts(&top_umi, &counts);
        assert_eq!(result, 15);
    }

    #[test]
    fn test_get_best_phred_single_group() {
        let mut clusters: Vec<Vec<Record>> = Vec::new();
        let mut reads = Vec::new();
        let mut record1 = Record::new();
        let mut record2 = Record::new();
        // Simulate BAM records with some sequences and quality scores
        record1.set(b"read1", None, b"ATCG", b"####");
        record2.set(b"read2", None, b"ATAG", b"##!#");
        reads.push(record1);
        reads.push(record2);
        clusters.push(reads);
        let result = get_best_phred(clusters);
        assert_eq!(result.seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_get_best_phred_multiple_groups() {
        let mut clusters: Vec<Vec<Record>> = Vec::new();
        let mut reads1 = Vec::new();
        let mut reads2 = Vec::new();
        let mut record1 = Record::new();
        let mut record2 = Record::new();
        let mut record3 = Record::new();
        let mut record4 = Record::new();
        // Simulate BAM records with some sequences and quality scores
        record1.set(b"read1", None, b"ATCG", b"####");
        record2.set(b"read2", None, b"ATCG", b"!!!!");
        record3.set(b"read3", None, b"ATTG", b"!!!!");
        record4.set(b"read4", None, b"ATGG", b"!!!!");
        reads1.push(record1);
        reads1.push(record2);
        reads2.push(record3);
        reads2.push(record4);
        clusters.push(reads1);
        clusters.push(reads2);
        let result = get_best_phred(clusters);
        assert_eq!(result.seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_correction_tied_for_majority() {
        let mut clusters: Vec<Vec<Record>> = Vec::new();
        let mut reads1 = Vec::new();
        let mut reads2 = Vec::new();
        let mut record1 = Record::new();
        let mut record2 = Record::new();
        // Simulate BAM records with some sequences and quality scores
        record1.set(b"read1", None, b"ATCG", b"####");
        record2.set(b"read2", None, b"ATAG", b"##!#");
        reads1.push(record1);
        reads2.push(record2);
        clusters.push(reads1);
        clusters.push(reads2);
        let result = get_best_phred(clusters);
        assert_eq!(result.seq().as_bytes(), b"ATCG");
    }
}
