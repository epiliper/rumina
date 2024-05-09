extern crate bam;
use crate::bottomhash::ReadsAndCount;
use crate::IndexMap;
use bam::Record;
use rand::Rng;

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper {
    // pub num: i32,
}

impl Grouper {
    // for umis successfully grouped (passing directional filtering)
    // remove the reads associated with each UMI from the bundle
    // tag them
    // push them to a list of tagged Records awaiting writing to an output bamfile
    pub fn tag_groups(
        &mut self,
        final_umis: Vec<Vec<&String>>,
        umis_records: &mut IndexMap<String, ReadsAndCount>,
    ) -> Vec<Record> {
        // for each UMI within a group, assign the same tag
        let mut rng = rand::thread_rng();
        let mut output_list: Vec<Record> = Vec::with_capacity(1_000_000);
        for top_umi in final_umis {
            // let ug_tag = self.num;
            let ug_tag = rng.gen_range(0..20_000_000);
            for group in top_umi {
                umis_records
                    .swap_remove(group)
                    .unwrap()
                    .reads
                    .drain(0..)
                    .for_each(|mut x| {
                        x.tags_mut().push_num(b"UG", ug_tag);
                        x.tags_mut().push_string(b"BX", group.as_bytes());
                        output_list.push(x);
                    })
            }
        }
        return output_list;
    }

    // driver function of grouping
    // recieves outout of main_grouper()
    pub fn tag_records(
        &mut self,
        grouping_output: Option<Vec<Vec<&String>>>,
        mut umis_records: IndexMap<String, ReadsAndCount>,
    ) -> Option<Vec<Record>> {
        match grouping_output {
            Some(groups) => {
                let reads = self.tag_groups(groups, &mut umis_records);
                return Some(reads);
            }
            None => return None,
        }
    }
}
