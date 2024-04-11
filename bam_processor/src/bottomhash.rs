extern crate bam;
extern crate rust_htslib;
use bam::Record;
use std::collections::HashMap;

pub struct BottomHashMap {
    pub bottom_dict: PositionKey,
}

impl BottomHashMap {
    pub fn update_dict(&mut self, position: &i32, key: i32, umi: &String, read: &Record) {
        self.bottom_dict
            .entry((*position).into())
            .or_default()
            .entry((key).into())
            .or_default()
            .entry(umi.into())
            .or_insert(ReadsAndCount {
                reads: Vec::new(),
                count: 0,
            })
            .up(read)
    }

    pub fn iter(&self) -> BundleIterator {
        BundleIterator {
            bottomhash: &self.bottom_dict,
            index: 0,
        }
    }
}

// this is the thing we iterate through
pub struct BundleIterator<'a> {
    bottomhash: &'a PositionKey,
    index: i32,
}

// this is what we return per iteration of said thing
impl<'a> Iterator for BundleIterator<'a> {
    type Item = (&'a UMIReads, i32);

    fn next(&mut self) -> Option<Self::Item> {
        while self.index < self.bottomhash.len().try_into().unwrap() {
            for p in self.bottomhash.keys() {
                for k in self.bottomhash[p].keys() {
                    return Some((&self.bottomhash[&p][&k], k.clone()));
                }
                self.index += 1;
            }
        }
        return None;
    }
}

type PositionKey = HashMap<i32, KeyUMI>; //every position has a key
type KeyUMI = HashMap<i32, UMIReads>; // every key has a UMI
pub type UMIReads = HashMap<String, ReadsAndCount>; // every UMI has a set of reads

#[derive(Debug)]
pub struct ReadsAndCount {
    pub reads: Vec<Record>,
    pub count: i32,
}

impl ReadsAndCount {
    fn up(&mut self, read: &Record) {
        self.count += 1;
        self.reads.push(read.clone());
    }
}
