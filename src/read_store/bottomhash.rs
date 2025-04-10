use crate::read_store::read_store::*;
use rust_htslib::bam::Record;

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
    pub read_dict: PositionKey,
}

impl BottomHashMap {
    pub fn update_dict(&mut self, position: i64, key: u64, umi: String, read: Record) {
        self.read_dict
            .entry(position)
            .or_default()
            .entry(key)
            .or_default()
            .entry(umi)
            .or_insert_with(|| ReadsAndCount {
                reads: Vec::new(),
                count: 0,
            })
            .up(read)
    }
}
