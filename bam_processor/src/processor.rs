extern crate bam;
extern crate rust_htslib;
use std::collections::HashMap;

use crate::bottomhash::UMIReads;

pub struct Processor {
    pub groups: UMIReads,
}

impl Processor {
    // create a list of indexes that split UMIs into substrings
    pub fn get_substr_slices(umi_length: i32, threshold: i32) -> Option<Vec<i32>> {
        let dividend = umi_length / threshold;
        let remainder = umi_length % threshold;
        let sub_sizes = vec![
            (dividend + 1) * remainder,
            (dividend * (threshold - remainder)),
        ];
        let mut offset = 0;
        let mut slices = Vec::new();

        if sub_sizes.is_empty() {
            return None;
        } else {
            for s in sub_sizes {
                slices.push(offset + s);
                offset += 1;
            }
            return Some(slices);
        }
    }

    pub fn build_substr_idx(
        mut umis: Vec<String>,
        umi_length: i32,
        threshold: i32,
    ) -> HashMap<i32, HashMap<String, Vec<String>>> {
        // use indices to split substrings, return a dict comparing all UMIs with similar codes
        let mut substr_idx: HashMap<i32, HashMap<String, Vec<String>>> = HashMap::new();
        let slices = Processor::get_substr_slices(umi_length, threshold + 1).unwrap();

        for idx in slices {
            for u in &mut umis {
                let u_sub = &u.split_off(idx.try_into().unwrap());
                substr_idx
                    .entry(idx)
                    .or_default()
                    .entry(u_sub.to_string())
                    .or_default()
                    .push(u.to_string());
            }
        }

        return substr_idx;
    }

    // pub fn iter_nearest_neighbors(umis: Vec<String>, substr_idx:HashMap<i32, HashMap<String, Vec<String>>>) {
    //     for (i, umi) in umis.iter().enumerate() {
    //         let neighbors: Vec<String> = Vec::new();

    //         for (idx, substring) in &substr_idx {

    //             u_sub = [umi.split_off(idx)z]
    //         }

    //     }

    // }
}
