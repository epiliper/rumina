extern crate bam;
extern crate rust_htslib;
use bam::Record;
use indexmap::IndexMap;

/* When main function executes, this struct is populated with
* all information necessary for grouping/deduplicating.
*/

/* This struct is a 4-layer deep nested hashmap.
* This is how it's organized, from top to bottom layer:
*
* positition along reference genome: {
*       key with optional metadata: {
*               all the umis at said position with said metadata {
*                   number of reads, and the Records themselves
*               }
*               }
* }
*/
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
}

type PositionKey = IndexMap<i32, KeyUMI>; //every position has a key
type KeyUMI = IndexMap<i32, UMIReads>; // every key has a UMI
pub type UMIReads = IndexMap<String, ReadsAndCount>; // every UMI has a set of reads

#[derive(Debug)]
pub struct ReadsAndCount {
    pub reads: Vec<Record>,
    pub count: i32,
}

// when a read is added to ReadsAndCount, increase the read count at x umi
impl ReadsAndCount {
    fn up(&mut self, read: &Record) {
        self.count += 1;
        self.reads.push(read.clone());
    }
}
