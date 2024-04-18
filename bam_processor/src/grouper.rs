extern crate bam;
use indexmap::IndexMap;
use crate::processor::Processor;

// this struct holds methods to
// 1. modify the records within the bottomhash by lookup
// 2. for every bundle, write the UG-tagged reads to output bam
pub struct Grouper{
}

impl Grouper {
    pub fn tag_by_group(
        adj_list: Vec<Vec<String>>,
        processor: Processor,
    ) {
        todo!();
    }
} 
