use crate::bottomhash::ReadsAndCount;
use crate::group_report::GroupReport;
use crate::grouper::GroupIterator;
use crate::read_picker::{correct_errors, get_counts, push_all_reads};
use crate::utils::{get_umi, get_umi_static};
use crate::IndexMap;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;

use core::str;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::collections::HashSet;

const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ\
                            abcdefghijklmnopqrstuvwxyz\
                            0123456789)(*&^%$#@!~";

const UMI_TAG_LEN: usize = 8;

// this struct serves to
// 1. for a given UMI group, pull all associated reads and
// 2. deduplicate by sequence majority or
// 3. output all reads in group
//
// remaining reads will be assigned a group-specific "UG" tag.
pub struct GroupHandler<'a> {
    pub seed: u64,
    pub group_only: bool,
    pub singletons: bool,
    pub separator: &'a String,
}

pub fn generate_tag(
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
        used_tags.insert(ug_tag);

        ug_tag
    } else {
        generate_tag(rng, used_tags)
    }
}

impl<'a> GroupHandler<'a> {
    // remove the reads associated with each UMI from the bundle
    // deduplicate and tag, or just tag them
    // push them to a list of tagged records awaiting writing to an output bamfile

    // driver function of grouping
    pub fn tag_records(
        &mut self,
        grouping_output: (HashMap<&str, i32>, Option<GroupIterator>),
        mut umis_records: IndexMap<String, ReadsAndCount>,
    ) -> Option<(Option<GroupReport>, Vec<Record>)> {
        match grouping_output.1 {
            Some(groups) => {
                let reads = self.tag_groups(groups, &mut umis_records, grouping_output.0);
                Some(reads)
            }
            None => None,
        }
    }
    pub fn tag_groups(
        &mut self,
        // mut final_umis: Vec<Vec<&str>>,
        final_umis: GroupIterator,
        umis_records: &mut IndexMap<String, ReadsAndCount>,
        counts: HashMap<&str, i32>,
    ) -> (Option<GroupReport>, Vec<Record>) {
        // for each UMI within a group, assign the same tag

        let mut rng = StdRng::seed_from_u64(self.seed);
        let mut output_list: Vec<Record> = Vec::new();
        let mut used_tags: HashSet<[u8; UMI_TAG_LEN]> = HashSet::with_capacity(output_list.len());

        // either group reads, or group and deduplicate
        let read_processor = match self.group_only {
            true => push_all_reads,
            false => correct_errors,
        };

        // to report min and max observed reads per group
        let mut group_report = GroupReport::new();

        let read_count_thres = match self.singletons {
            true => 0,
            false => 3, // groups should have 3+ reads for reliable majority rule
        };

        for top_umi in final_umis {
            let num_reads_in_group = get_counts(&top_umi, &counts);
            group_report.num_groups += 1;
            if num_reads_in_group >= read_count_thres {
                let ug_tag = generate_tag(&mut rng, &mut used_tags);

                // check if number of reads per group is new minimum or maximum
                if num_reads_in_group < group_report.min_reads_per_group {
                    group_report.min_reads_per_group = num_reads_in_group;
                    group_report.min_reads_group = ug_tag;
                }

                if num_reads_in_group > group_report.max_reads_per_group {
                    group_report.max_reads_per_group = num_reads_in_group;
                    group_report.max_reads_group = ug_tag;
                }

                // since the group has enough reads to be used, count it in the report
                group_report.num_passing_groups += 1;

                let mut cluster_list: Vec<ReadsAndCount> = Vec::new();

                // get all the reads across all the umis in the group
                for group in &top_umi {
                    cluster_list.push(umis_records.swap_remove(*group).unwrap());
                }

                // tag final reads and send for writing to output bam
                let mut to_write = read_processor(&mut cluster_list);

                to_write.iter_mut().for_each(|read| {
                    let read_umi = get_umi_static(get_umi(read, self.separator));

                    // add group tag
                    read.push_aux(b"UG", Aux::String(str::from_utf8(&ug_tag).unwrap()))
                        .unwrap();

                    read.push_aux(
                        b"BX",
                        Aux::String(str::from_utf8(read_umi.as_slice()).unwrap()),
                    )
                    .unwrap();

                    group_report.num_reads_output_file += 1;
                });

                output_list.extend(to_write);
            }
        }

        if group_report.is_blank() {
            (None, output_list)
        } else {
            (Some(group_report), output_list)
        }
    }
}
