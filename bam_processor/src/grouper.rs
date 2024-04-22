extern crate bam;
use crate::bottomhash::BottomHashMap;
use rand::{thread_rng, Rng};

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper {
}

impl Grouper {
    pub fn modify_records(
        &self,
        final_umis: Vec<Vec<String>>,
        mut bottomhash: BottomHashMap,
        position: i32,
    ) {
        let mut rng = thread_rng();

        for top_umi in final_umis {
            for group in top_umi {

                 bottomhash.bottom_dict[&position][0].entry(group).and_modify(
                |e| {
                        e.reads.iter_mut().for_each(
                            |x| x.tags_mut().push_num(b"UG", rng.gen_range(1_000_000..1_999_999))
                        );
                        });
                    }
            }
        }
    }
