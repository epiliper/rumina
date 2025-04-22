use crate::readkey::ReadKey;
use anyhow::{Context, Error};
use bio::io::fastq;
use core::str;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux};
use smol_str::SmolStr;

/// The record interface serves to allow using FASTQ and BAM records with the clustering and
/// deduplication portions of the pipeline. The [FastqRecord] and [BamRecord] structs are very thin
/// wrappers around their original types.
pub trait Record {
    fn seq(&self) -> String;
    fn seq_str(&self) -> &[u8];
    fn qual(&self) -> &[u8];
    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error>;
    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey);
    fn mark_group(&mut self, tag: &[u8]);
}

/// A wrapper around [rust_htslib::bam::Record]
pub type BamRecord = bam::Record;

impl Record for BamRecord {
    fn seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq().encoded.to_vec()) }
    }

    fn seq_str(&self) -> &[u8] {
        self.seq().encoded
    }

    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error> {
        // unsafe {
        //     Ok(std::str::from_utf8_unchecked(self.qname())
        //         .rsplit_once(separator)
        //         .context("failed to get UMI from read QNAME. Check --separator. Exiting.")?
        //         .1
        //         .to_string())
        // }
        //
        unsafe {
            let s = std::str::from_utf8_unchecked(self.qname())
                .rsplit_once(separator)
                .with_context(|| {
                    format!(
                        "failed to get UMI with separator {}. QNAME in question:\n{}",
                        separator,
                        std::str::from_utf8(self.qname()).unwrap()
                    )
                })?
                .1;

            Ok(SmolStr::from(s))
        }
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey) {
        let mut pos;
        let key: ReadKey;

        if self.is_reverse() {
            // set end pos as start to group with forward-selfs covering same region
            pos = self.reference_end();
            pos += self.cigar().trailing_softclips(); // pad with right-side soft clip
            key = ReadKey {
                length: self.seq_len() * group_by_length as usize,
                reverse: true,
                chr: self.tid() as usize,
            };
            (pos, key)
        } else {
            pos = self.reference_start();
            pos -= self.cigar().leading_softclips(); // pad with left-side soft clip

            key = ReadKey {
                length: self.seq_len() * group_by_length as usize,
                reverse: false,
                chr: self.tid() as usize,
            };
            (pos, key)
        }
    }

    fn mark_group(&mut self, tag: &[u8]) {
        self.push_aux(b"BX", Aux::String(str::from_utf8(tag).unwrap()))
            .unwrap();
    }
}

/// A wrapper around [bio::io::fastq::Record]
// pub type FastqRecord = fastq::Record;
pub type FastqRecord = seq_io::fastq::OwnedRecord;

impl Record for FastqRecord {
    fn seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq.clone()) }
    }

    fn seq_str(&self) -> &[u8] {
        self.seq.as_slice()
    }

    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error> {
        unsafe {
            let s = std::str::from_utf8_unchecked(&self.head)
                .rsplit_once(separator)
                .with_context(|| {
                    format!(
                        "failed to get UMI with separator {}. QNAME in question:\n{}",
                        separator,
                        self.seq()
                    )
                })?
                .1;

            Ok(SmolStr::from(s))
        }
    }

    fn qual(&self) -> &[u8] {
        self.qual.as_slice()
    }

    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey) {
        let pos = 1;
        let key = ReadKey {
            length: self.seq().len() & group_by_length as usize,
            reverse: false,
            chr: 1,
        };

        (pos, key)
    }

    fn mark_group(&mut self, _tag: &[u8]) {}
}
