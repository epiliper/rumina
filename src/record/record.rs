use crate::readkey::ReadKey;
use anyhow::{Context, Error};
use core::str;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux};
use seq_io::fastq::Record as _;
use smol_str::SmolStr;

pub fn extract_umi_from_header<'a>(header: &'a str, separator: &String) -> Result<&'a str, Error> {
    let (_rest, past_sep) = header.rsplit_once(separator).with_context(|| {
        format!(
            "failed to get UMI with separator '{}'. Header in question:\n{}",
            separator, header
        )
    })?;

    // check if we have r1/r2 extensions, e.g:
    // <REST_OF_HEADER>_<UMI>/1
    // if we don't remove that, we'll over-stratify read groups.
    if let Some((umi, _mate_info)) = past_sep.rsplit_once("/") {
        Ok(umi)
    } else {
        Ok(past_sep)
    }
}

pub fn reverse_complement(s: &str) -> SmolStr {
    let mut out: Vec<u8> = Vec::with_capacity(s.len());
    for c in s.bytes().rev() {
        let r = match c.to_ascii_uppercase() {
            b'A' => b'T',
            b'G' => b'C',
            b'C' => b'G',
            b'T' => b'A',
            b'N' => b'N',
            _ => panic!(),
        };

        out.push(r);
    }

    unsafe {
        let out = std::str::from_utf8_unchecked(&out);
        SmolStr::new(&out)
    }
}

/// The record interface serves to allow using FASTQ and BAM records with the clustering and
/// deduplication portions of the pipeline. The [FastqRecord] and [BamRecord] structs are very thin
/// wrappers around their original types.
pub trait Record {
    fn _seq(&self) -> String;
    fn seq_str(&self) -> &[u8];
    fn qual(&self) -> &[u8];
    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error>;
    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey);
    fn mark_group(&mut self, tag: &[u8]);
}

/// A wrapper around [rust_htslib::bam::Record]
pub type BamRecord = bam::Record;

impl Record for BamRecord {
    fn _seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq().encoded.to_vec()) }
    }

    fn seq_str(&self) -> &[u8] {
        self.seq().encoded
    }

    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error> {
        unsafe {
            let s = std::str::from_utf8_unchecked(self.qname());
            Ok(SmolStr::from(extract_umi_from_header(s, separator)?))
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
    fn _seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq.clone()) }
    }

    fn seq_str(&self) -> &[u8] {
        self.seq.as_slice()
    }

    fn get_umi(&self, separator: &String) -> Result<SmolStr, Error> {
        Ok(SmolStr::from(extract_umi_from_header(
            self.id()?,
            separator,
        )?))
    }

    fn qual(&self) -> &[u8] {
        self.qual.as_slice()
    }

    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey) {
        let pos = 1;
        let key = ReadKey {
            length: self.seq_str().len() & group_by_length as usize,
            reverse: false,
            chr: 1,
        };

        (pos, key)
    }

    fn mark_group(&mut self, _tag: &[u8]) {}
}

#[test]
fn test_rev1() {
    let o = "GTCTATATA";
    assert_eq!(reverse_complement(o), "TATATAGAC");
    assert_eq!(reverse_complement(&reverse_complement(o)), o);
}
