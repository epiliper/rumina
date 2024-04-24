extern crate bam;
use crate::IndexMap;
use bam::Record;
use crate::bottomhash::ReadsAndCount;
use crate::processor::GroupsAndSingletons;
use rand::{thread_rng, Rng};

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper {}

impl Grouper {
    pub fn tag_groups(
        &self,
        final_umis: Vec<Vec<&String>>,
        // bottomhash: & mut BottomHashMap,
        umis_records: & mut IndexMap<String, ReadsAndCount>,
        output_list: & mut Vec<Record>,
    ) {
        let mut rng = thread_rng();

        // for each UMI within a group, assign the same tag
        for top_umi in final_umis {
            for group in top_umi {
                let ug_tag = rng.gen_range(1_000_000..10_999_999);

                   umis_records.swap_remove(group).unwrap().reads.drain(0..)
                .for_each(
                        |mut x| {
                            x.tags_mut().push_num(b"UG", ug_tag);
                            output_list.push(x);
                        }
                    )
                }
                // umis_records.entry(group.to_string())
                //     .and_modify(|e| {
                //         e.reads
                //             .iter_mut()
                //             .for_each(|x| {
                //                 x.tags_mut().push_num(b"UG", ug_tag);
                //                 output_list.push(*x);
                //             }
                //             );
                //     });
            }
        }

    pub fn tag_singletons(
        &self,
        singletons: Vec<&String>,
        // bottomhash: & mut BottomHashMap,
        umis_records: & mut IndexMap<String, ReadsAndCount>,
        output_list: & mut Vec<Record>,
    ) {
        println!{"YYYYY"};
        let mut rng = thread_rng();
        for dud in singletons {
            println!{"{}", dud};
            let ug_tag = rng.gen_range(1_000_000..10_999_999);
            umis_records.swap_remove(dud).unwrap().reads.drain(0..)
                .for_each(
                        |mut x| {
                            x.tags_mut().push_num(b"UG", ug_tag);
                            output_list.push(x);
                        }
                    )
                }
            // umis_records.get_mut(dud).expect("UMI not found in sequence!")
            // .reads
            // .iter_mut()
            //     .for_each(|x| x.tags_mut().push_num(b"UG", ug_tag));
    }

    pub fn tag_records(
        &self,
        grouping_output: GroupsAndSingletons,
        mut umis_records: IndexMap<String, ReadsAndCount>,
        output_list: & mut Vec<Record>,
    ) {
        match grouping_output {
            (Some(groups), Some(singletons)) => {
                self.tag_groups(groups, &mut umis_records, output_list);
                self.tag_singletons(singletons, &mut umis_records, output_list);
            }

            (None, Some(singletons)) => {
                self.tag_singletons(singletons, &mut umis_records, output_list);
            }

            (Some(groups), None) => {
                self.tag_groups(groups, &mut umis_records, output_list);
            }

            (None, None) => {}
        }
    }

}
