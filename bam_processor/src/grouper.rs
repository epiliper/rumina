extern crate bam;
use crate::bottomhash::ReadsAndCount;
use crate::IndexMap;
use bam::Record;
use rand::rngs::ThreadRng;
use rand::Rng;
use std::collections::HashMap;
use std::collections::HashSet;
use crate::dedup_correct::{get_counts, correct_errors};


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


impl Grouper {
    // remove the reads associated with each UMI from the bundle
    // deduplicate and error correct them
    // tag filtered reads
    // push them to a list of tagged records awaiting writing to an output bamfile
    pub fn tag_groups(
        &mut self,
        final_umis: Vec<Vec<&String>>,
        umis_records: &mut IndexMap<String, ReadsAndCount>,
        counts: HashMap<&String, i32>,
    ) -> (Option<IndexMap<i64, i64>>, Vec<Record>) {
        // for each UMI within a group, assign the same tag

        let mut rng = rand::thread_rng();
        let mut output_list: Vec<Record> = Vec::with_capacity(1_000_000);
        let mut used_tags: HashSet<[u8; UMI_TAG_LEN]> = HashSet::with_capacity(output_list.len());

        let mut first = true;

        // to report min and max observed reads per group
        let mut min_max_report: IndexMap<i64, i64> = IndexMap::new();

        for top_umi in final_umis {
            let num_reads_in_group = get_counts(&top_umi, &counts);
            if num_reads_in_group >= 3 {

                if first {
                    // initialize min and max for comparison
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

                // get all the reads across all the umis in the group
                for group in &top_umi {
                    cluster_list.push(umis_records.swap_remove(*group).unwrap());
                }

                // this is the highest qual read from the read group with the highest average phred score
                let mut final_record = correct_errors(&mut cluster_list);

                final_record.tags_mut().push_string(b"UG", &ug_tag);
                final_record.tags_mut().push_string(b"BX", &top_umi.iter().next().unwrap().as_bytes());
                output_list.push(final_record);
                }
            }

        if min_max_report.is_empty() {
            return (None, output_list);
        }
        return (Some(min_max_report), output_list);
    }
    

    // driver function of grouping
    // recieves outout of main_grouper()
    pub fn tag_records(
        &mut self,
        grouping_output: (HashMap<&String, i32>, Option<Vec<Vec<&String>>>),
        mut umis_records: IndexMap<String, ReadsAndCount>,
    ) -> Option<(Option<IndexMap<i64, i64>>, Vec<Record>)> {

        match grouping_output.1 {
            Some(groups) => {
                let reads = self.tag_groups(groups, &mut umis_records, grouping_output.0);
                return Some(reads);
            }
            None => return None,
        }
    }
}
