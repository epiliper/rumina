use crate::group_report::GroupReport;
use crate::processor::UmiHistogram;
use crate::read_picker::{correct_errors, get_counts, push_all_reads};
use crate::read_store::{ReadStore, SeqMap};
use crate::record::Record;
use crate::IndexMap;
use indexmap::IndexSet;

use anyhow::{Context, Error, Result};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
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
    pub min_depth: usize,
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
    pub fn tag_records<T: Record>(
        &mut self,
        // mut final_umis: Vec<Vec<&str>>,
        final_umis: impl Iterator<Item = IndexSet<smol_str::SmolStr>>,
        umis_records: &mut IndexMap<smol_str::SmolStr, (i32, SeqMap<T>)>,
        counts: UmiHistogram,
    ) -> Result<(Option<GroupReport>, Vec<T>), Error> {
        // for each UMI within a group, assign the same tag

        let mut rng = StdRng::seed_from_u64(self.seed);
        let mut output_list: Vec<T> = Vec::new();
        let mut used_tags: HashSet<[u8; UMI_TAG_LEN]> = HashSet::with_capacity(output_list.len());

        // either group reads, or group and deduplicate
        let read_processor = match self.group_only {
            true => push_all_reads,
            false => correct_errors,
        };

        // to report min and max observed reads per group
        let mut group_report = GroupReport::new();

        // let read_count_thres = match self.min_depth {
        //     true => 0,
        //     false => 3, // groups should have 3+ reads for reliable majority rule
        // };
        //
        let read_count_thres = self.min_depth;

        for top_umi in final_umis {
            let num_reads_in_group = get_counts(&top_umi, &counts);
            group_report.num_groups += 1;
            if num_reads_in_group >= read_count_thres as i64 {
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

                let mut top_umi = top_umi.iter();
                let (_, mut seq_map) = umis_records
                    .swap_remove(top_umi.next().unwrap())
                    .context("FATAL: Attempted to remove UMI that was already sent for output")?;

                // get all the reads across all the umis in the group
                for group in top_umi {
                    seq_map.combine(umis_records.swap_remove(group).unwrap().1, self.group_only);
                }

                // tag final reads and send for writing to output bam
                let mut to_write = read_processor(&mut seq_map);

                // TODO: FIX THIS FOR BOTH BAM AND FASTQ RECORDS
                to_write.iter_mut().for_each(|read| {
                    // TODO; avoid clone
                    if let Some(umi) = read.get_umi(self.separator).ok() {
                        read.mark_group(umi.as_bytes());
                    } else {
                        println! {"Unable to retrieve UMI for output read! Skipping..."}
                    }

                    // add group tag
                    // read.push_aux(b"UG", Aux::String(str::from_utf8(&ug_tag).unwrap()))
                    //     .unwrap();

                    group_report.num_reads_output_file += 1;
                });

                output_list.extend(to_write);
            }
        }

        if group_report.is_blank() {
            Ok((None, output_list))
        } else {
            Ok((Some(group_report), output_list))
        }
    }
}
