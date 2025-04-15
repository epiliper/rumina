use indexmap::IndexMap;
// use rust_htslib::bam::Record;
use crate::record::record::Record;

pub type PositionKey<T> = IndexMap<i64, KeyUMI<T>>; //every position has a key
type KeyUMI<T> = IndexMap<u64, UMIReads<T>>; // every key has a UMI
pub type UMIReads<T> = IndexMap<String, ReadsAndCount<T>>; // every UMI has a set of reads

#[derive(Debug)]
pub struct ReadsAndCount<T>
where
    T: Record,
{
    pub reads: Vec<T>,
    pub count: i32,
}

// when a read is added to ReadsAndCount, increase the read count at x umi
impl<T: Record> ReadsAndCount<T> {
    pub fn up(&mut self, read: T) {
        self.count += 1;
        self.reads.push(read);
    }
}
