use crate::io::fastqio::FastqRecord;
use crate::readkey::ReadKey;
use anyhow::{Context, Error};
use core::str;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux};
use smol_str::SmolStr;

pub fn extract_umi_from_header<'a>(header: &'a str, separator: &str) -> Result<&'a str, Error> {
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
        SmolStr::new(out)
    }
}

/// The record interface serves to allow using FASTQ and BAM records with the clustering and
/// deduplication portions of the pipeline. The [FastqRecord] and [BamRecord] structs are very thin
/// wrappers around their original types.
pub trait SequenceRecord {
    fn _seq(&self) -> String;
    fn seq_str(&self) -> &[u8];
    fn qual(&self) -> &[u8];
    fn get_umi(&self, separator: &str) -> Result<SmolStr, Error>;
    fn get_pos_key(&self, group_by_length: bool) -> (i64, ReadKey);
    fn mark_group(&mut self, umi: &[u8], group_tag: &[u8]);
    fn qname(&self) -> &[u8];
}

/// A wrapper around [rust_htslib::bam::Record]
pub type BamRecord = bam::Record;

impl SequenceRecord for BamRecord {
    fn _seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq().encoded.to_vec()) }
    }

    fn qname(&self) -> &[u8] {
        self.qname()
    }

    fn seq_str(&self) -> &[u8] {
        self.seq().encoded
    }

    fn get_umi(&self, separator: &str) -> Result<SmolStr, Error> {
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
                // using seq_len_from_cigar to consider soft-clipped bases, essentially getting
                // the entire read length. Not going to consider hard-clipping for now.
                length: self.seq_len_from_cigar(false) * group_by_length as usize,
                reverse: true,
                chr: self.tid() as usize,
            };
            (pos, key)
        } else {
            pos = self.reference_start();
            pos -= self.cigar().leading_softclips(); // pad with left-side soft clip

            key = ReadKey {
                length: self.seq_len_from_cigar(false) * group_by_length as usize,
                reverse: false,
                chr: self.tid() as usize,
            };
            (pos, key)
        }
    }

    fn mark_group(&mut self, umi: &[u8], group_tag: &[u8]) {
        self.push_aux(b"BX", Aux::String(str::from_utf8(umi).unwrap()))
            .unwrap();
        self.push_aux(b"UG", Aux::String(str::from_utf8(group_tag).unwrap()))
            .unwrap();
    }
}

/// A wrapper around [bio::io::fastq::Record]
// pub type FastqRecord = fastq::Record;
impl SequenceRecord for FastqRecord {
    fn _seq(&self) -> String {
        unsafe { std::str::from_utf8_unchecked(self.seq()).to_string() }
    }

    fn seq_str(&self) -> &[u8] {
        self.seq()
    }

    fn get_umi(&self, separator: &str) -> Result<SmolStr, Error> {
        Ok(SmolStr::from(extract_umi_from_header(
            // self.id()?,
            self.id(),
            separator,
        )?))
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn qname(&self) -> &[u8] {
        self.id().as_bytes()
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

    fn mark_group(&mut self, _tag: &[u8], _group_tag: &[u8]) {}
}

#[test]
fn test_rev1() {
    let o = "GTCTATATA";
    assert_eq!(reverse_complement(o), "TATATAGAC");
    assert_eq!(reverse_complement(&reverse_complement(o)), o);
}
