use crate::read_store::read_store::*;
use crate::record::Record;
use smol_str::SmolStr;

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
pub struct BottomHashMap<T: Record> {
    pub read_dict: PositionKey<T>,
    pub read_count: u64,
    // pub update_method: fn(read: T),
}

impl<T: Record> BottomHashMap<T> {
    pub fn update_dict(
        &mut self,
        position: i64,
        key: u64,
        umi: SmolStr,
        read: T,
        retain_all: bool,
    ) {
        let (count, seq_map) = self
            .read_dict
            .entry(position)
            .or_default()
            .entry(key)
            .or_default()
            .entry(umi)
            .or_insert_with(|| (0, SeqMap::new()));

        *count += 1;
        self.read_count += seq_map.intake(read, retain_all) as u64;
    }

    pub fn shrink_to_fit(&mut self) {
        self.read_dict.iter_mut().for_each(|(_pos, key_map)| {
            key_map.shrink_to_fit();
            key_map.iter_mut().for_each(|(_key, seq_dict)| {
                seq_dict.shrink_to_fit();
            })
        })
    }
}
