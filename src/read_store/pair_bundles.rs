use crate::read_store::read_store::*;
use crate::record::BamRecord;
use indexmap::IndexMap;
use log::warn;
use rust_htslib::bam::{record::Aux, Record};

pub struct PairBundles {
    pub read_dict: IndexMap<String, ReadsAndCount<BamRecord>>,
}

const UMI_TAG: &[u8; 2] = b"BX";

impl PairBundles {
    pub fn update_dict(&mut self, read: Record) {
        let umi = if let Ok(Aux::String(bx_i)) = read.aux(UMI_TAG) {
            bx_i
        } else {
            warn!("Cannot find UMI for read: {:?}", read);
            "NULL"
        };

        self.read_dict
            .entry(umi.to_string())
            .or_insert_with(|| ReadsAndCount {
                reads: Vec::new(),
                count: 0,
            })
            .up(read)
    }
}
