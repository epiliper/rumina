use indexmap::IndexMap;
use rust_htslib::bam::Record;

pub type PositionKey = IndexMap<i64, KeyUMI>; //every position has a key
type KeyUMI = IndexMap<u64, UMIReads>; // every key has a UMI
pub type UMIReads = IndexMap<String, ReadsAndCount>; // every UMI has a set of reads

#[derive(Debug)]
pub struct ReadsAndCount {
    pub reads: Vec<Record>,
    pub count: i32,
}

// when a read is added to ReadsAndCount, increase the read count at x umi
impl ReadsAndCount {
    pub fn up(&mut self, read: Record) {
        self.count += 1;
        self.reads.push(read);
    }
}
