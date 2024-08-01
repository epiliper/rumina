extern crate bam;
use crate::bottomhash::ReadsAndCount;
use crate::read_io::GroupReport;
use crate::read_picker::{correct_errors, get_counts};
use crate::IndexMap;
use bam::Record;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::collections::HashSet;

const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ\
                            abcdefghijklmnopqrstuvwxyz\
                            0123456789)(*&^%$#@!~";

const UMI_TAG_LEN: usize = 8;

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Deduplicator {
    pub seed: u64,
}

pub fn generate_tag(
    // rng: &mut ThreadRng,
    rng: &mut StdRng,
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

impl Deduplicator {
    // remove the reads associated with each UMI from the bundle
    // deduplicate and error correct them
    // tag filtered reads
    // push them to a list of tagged records awaiting writing to an output bamfile

    // driver function of grouping
    pub fn tag_records(
        &mut self,
        grouping_output: (HashMap<&String, i32>, Option<Vec<Vec<&String>>>),
        mut umis_records: IndexMap<String, ReadsAndCount>,
    ) -> Option<(Option<GroupReport>, Vec<Record>)> {
        match grouping_output.1 {
            Some(groups) => {
                let reads = self.tag_groups(groups, &mut umis_records, grouping_output.0);
                return Some(reads);
            }
            None => return None,
        }
    }
    pub fn tag_groups(
        &mut self,
        mut final_umis: Vec<Vec<&String>>,
        umis_records: &mut IndexMap<String, ReadsAndCount>,
        counts: HashMap<&String, i32>,
    ) -> (Option<GroupReport>, Vec<Record>) {
        // for each UMI within a group, assign the same tag

        let mut rng = StdRng::seed_from_u64(self.seed);
        let mut output_list: Vec<Record> = Vec::with_capacity(1_000_000);
        let mut used_tags: HashSet<[u8; UMI_TAG_LEN]> = HashSet::with_capacity(output_list.len());

        let mut first = true;

        // to report min and max observed reads per group
        let mut group_report: GroupReport = Default::default();
        group_report.num_groups = 0;
        group_report.num_groups += final_umis.len() as i64;

        for top_umi in final_umis.drain(0..) {
            let num_reads_in_group = get_counts(&top_umi, &counts);
            if num_reads_in_group >= 3 {
                let ug_tag = generate_tag(&mut rng, &mut used_tags);

                if first {
                    // initialize min and max for comparison
                    first = false;
                    group_report.min_reads = i64::MAX;
                    group_report.min_reads_group = *b"NONENONE";

                    group_report.max_reads = 0;
                    group_report.max_reads_group = *b"NONENONE";
                }
                // check if number of reads per group is new minimum or maximum
                if num_reads_in_group < group_report.min_reads {
                    group_report.min_reads = num_reads_in_group;
                    group_report.min_reads_group = ug_tag;
                }

                if num_reads_in_group > group_report.max_reads {
                    group_report.max_reads = num_reads_in_group;
                    group_report.max_reads_group = ug_tag;
                }

                // since the group has enough reads to be used, count it in the report
                group_report.num_passing_groups += 1;

                let mut cluster_list: Vec<ReadsAndCount> = Vec::new();

                // get all the reads across all the umis in the group
                for group in &top_umi {
                    cluster_list.push(umis_records.swap_remove(*group).unwrap());
                }

                // this is the read from the majority sequence read group
                let mut final_record = correct_errors(&mut cluster_list);

                final_record.tags_mut().push_string(b"UG", &ug_tag);
                final_record
                    .tags_mut()
                    .push_string(b"BX", &top_umi.iter().next().unwrap().as_bytes());
                output_list.push(final_record);
            }
        }

        if first {
            return (None, output_list);
        }
        return (Some(group_report), output_list);
    }
}
