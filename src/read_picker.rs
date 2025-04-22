/* This module contains scripts for deduplication and error-correction
* of reads within a UMI group.
*
* correct_errors() is the driver function;
* 1. input all the reads within the UMI group
* 2. group the reads by sequence
* 3. select the group with the highest phred score across reads
* 4. Output one read from the group
*/

use crate::read_store::read_store::SeqMap;
use crate::record::Record;
use indexmap::IndexSet;
use std::collections::HashMap;

pub fn correct_errors<T: Record>(clusters: &mut SeqMap<T>) -> Vec<T> {
    assert!(!clusters.is_empty());
    // sort, in descending order, by:
    // 1. sequence count, to get majority sequence
    // 2. summed quality score, to break potential ties from 1).
    clusters.sort_by(|_seq1, seq_entry1, _seq2, seq_entry2| {
        seq_entry2
            .count
            .cmp(&seq_entry1.count)
            .then(seq_entry2.qual_sum.cmp(&seq_entry1.qual_sum))
    });

    // return the first read of the sequence cluster at the top of the sort
    let (_, mut seq_entry) = clusters.swap_remove_index(0).unwrap();
    vec![seq_entry.reads.swap_remove(0)]
}

// used with the --group_only arg to return all reads within a group with a group tag
pub fn push_all_reads<T: Record>(clusters: &mut SeqMap<T>) -> Vec<T> {
    let mut reads_to_write: Vec<T> = Vec::with_capacity(clusters.len());

    clusters
        .drain(..)
        .for_each(|(_seq, seq_entry)| reads_to_write.extend(seq_entry.reads));

    reads_to_write
}

// get the number of reads across all UMIs within a group
// this is useful for setting a threshold for reads observed per UMI group
pub fn get_counts(top_umi: &IndexSet<smol_str::SmolStr>, counts: &HashMap<&str, i32>) -> i64 {
    let mut read_count = 0;
    for umi in top_umi {
        read_count += counts.get(umi.as_str()).unwrap();
    }
    read_count as i64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_store::read_store::{ReadStore, SeqMap};
    use rust_htslib::bam::Record;
    use std::collections::HashMap;

    #[test]
    fn test_correct_errors_single_read() {
        let mut record = Record::new();
        record.set(b"read1", None, b"ATCG", b"####");

        let mut cluster = SeqMap::new();
        cluster.intake(record, true);

        let result = correct_errors(&mut cluster);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_correct_errors_multiple_reads() {
        let mut cluster: SeqMap<Record> = SeqMap::new();
        // let mut reads1 = Vec::new();
        // let mut reads2 = Vec::new();
        let mut record1 = Record::new();
        let mut record2 = Record::new();
        let mut record3 = Record::new();
        // Simulate BAM records with some sequences and quality scores
        record1.set(b"read1", None, b"ATCG", b"####");
        record2.set(b"read2", None, b"ATCG", b"####");
        record3.set(b"read3", None, b"ATGG", b"####");

        cluster.intake(record1, true);
        cluster.intake(record2, true);
        cluster.intake(record3, true);

        let result = correct_errors(&mut cluster);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].seq().as_bytes(), b"ATCG");
    }

    #[test]
    fn test_get_counts() {
        let top_umi = IndexSet::from([
            smol_str::SmolStr::new("UMI1"),
            smol_str::SmolStr::new("UMI2"),
        ]);
        let counts = HashMap::from([("UMI1", 10), ("UMI2", 5)]);
        let result = get_counts(&top_umi, &counts);
        assert_eq!(result, 15);
    }
}
