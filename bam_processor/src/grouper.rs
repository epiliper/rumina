extern crate bam;
use crate::bottomhash::BottomHashMap;
use crate::processor::GroupsAndSingletons;
use rand::{thread_rng, Rng};

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper {}

impl Grouper {
    pub fn tag_groups(
        &self,
        final_umis: Vec<Vec<String>>,
        bottomhash: & mut BottomHashMap,
        position: i32,
    ) {
        let mut rng = thread_rng();

        // for each UMI within a group, assign the same tag
        for top_umi in final_umis {
            for group in top_umi {
                let ug_tag = rng.gen_range(1_000_000..10_999_999);
                bottomhash.bottom_dict[&position][0]
                    .entry(group)
                    .and_modify(|e| {
                        e.reads
                            .iter_mut()
                            .for_each(|x| x.tags_mut().push_num(b"UG", ug_tag));
                    });
            }
        }
    }

    pub fn tag_singletons(
        &self,
        singletons: Vec<&String>,
        bottomhash: & mut BottomHashMap,
        position: i32,
    ) {
        let mut rng = thread_rng();
        for dud in singletons {
            let ug_tag = rng.gen_range(1_000_000..10_999_999);
            // bottomhash.bottom_dict[&position][0]
            //     .entry(dud)
            //     .and_modify(|e| {
            //         e.reads
            //             .iter_mut()
            //             .for_each(|x| x.tags_mut().push_num(b"UG", ug_tag));
            //     });
            bottomhash.bottom_dict[&position][0].get_mut(dud).expect("UMI not found in sequence!")
            .reads
            .iter_mut()
                .for_each(|x| x.tags_mut().push_num(b"UG", ug_tag));
        }
    }

    pub fn tag_records(
        &self,
        grouping_output: GroupsAndSingletons,
        bottomhash: & mut BottomHashMap,
        position: i32,
    ) {
        match grouping_output {
            (Some(groups), Some(singletons)) => {
                self.tag_groups(groups, bottomhash, position);
                self.tag_singletons(singletons, bottomhash, position);
            }

            (None, Some(singletons)) => {
                self.tag_singletons(singletons, bottomhash, position);
            }

            (Some(groups), None) => {
                self.tag_groups(groups, bottomhash, position);
            }

            (None, None) => {}
        }
    }
}
