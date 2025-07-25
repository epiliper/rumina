use crate::record::record::BamRecord;
use indexmap::IndexMap;
use log::{debug, warn};
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use std::str;
use std::sync::Arc;

use crate::merge_report::*;
use crate::read_store::read_store::ReadsAndCount;
use crate::realign::{align_to_ref, ReMapper};

pub fn handle_dupes(
    umis_reads: &mut IndexMap<String, ReadsAndCount<BamRecord>>,
    mapper: ReMapper,
    ref_fasta: &Vec<u8>,
    min_overlap_bp: usize,
    sender: crossbeam::channel::Sender<Option<BamRecord>>,
) -> Vec<MergeResult> {
    let results: Arc<Mutex<Vec<MergeResult>>> = Arc::new(Mutex::new(Vec::new()));

    // for (_umi, mut reads) in umis_reads.drain() {
    umis_reads
        .par_drain(..)
        .for_each(|(_umi, reads_and_count)| {
            let mut reads = reads_and_count.reads;
            match reads.len() {
                1 => {
                    reads
                        .drain(..)
                        .for_each(|read| sender.send(Some(read)).unwrap());
                }
                _ => {
                    debug!("MERGING: number of reads with same umi: {}", reads.len());
                    // let mut outreads: Vec<Record> = Vec::with_capacity(100);
                    let mut merge_results: Vec<MergeResult> = Vec::with_capacity(50);

                    // sort the read list by reads likely to overlap
                    reads.sort_by(|ra, rb| {
                        ra.tid()
                            .cmp(&rb.tid())
                            .then_with(|| ra.qname().cmp(&rb.qname()))
                            .then_with(|| ra.is_reverse().cmp(&!rb.is_reverse()))
                            .then_with(|| ra.pos().cmp(&rb.pos()))
                    });

                    while !reads.is_empty() {
                        let read = reads.remove(0);

                        let result = find_merges(&read, &mut reads, min_overlap_bp);

                        match result {
                            MergeResult::Discordant(_) => {
                                merge_results.push(result);
                            }

                            MergeResult::NoMerge(_) => {
                                sender.send(Some(read)).unwrap();
                                merge_results.push(result);
                            }
                            MergeResult::Merge(merged_bases) => {
                                let (_start_pos, merged_seq) =
                                    construct_sequence(merged_bases.unwrap());
                                let merged_read = construct_read(
                                    &read,
                                    merged_seq,
                                    &mut mapper.clone(),
                                    &ref_fasta,
                                );
                                sender.send(Some(merged_read)).unwrap();
                                merge_results.push(MergeResult::Merge(None));
                            }
                        }
                    }
                    results.lock().extend(merge_results);
                }
            }
        });

    Arc::try_unwrap(results).unwrap().into_inner()
}

pub fn is_opp_orientation(read_a: &BamRecord, read_b: &BamRecord) -> bool {
    read_a.is_reverse() | read_b.is_reverse()
}

pub fn is_overlap(read_a: &BamRecord, read_b: &BamRecord) -> bool {
    let (ras, rae) = (read_a.reference_start(), read_a.reference_end());
    let (rbs, rbe) = (read_b.reference_start(), read_b.reference_end());

    if rae == rbe && ras == rbs {
        return true;
    }

    ras < rbs && rae >= rbs
}

// for groups of >2 reads, find every overlapping f/r read pair, attempt merge
pub fn find_merges(
    read: &BamRecord,
    reads: &mut Vec<BamRecord>,
    min_overlap_bp: usize,
) -> MergeResult {
    for (i, other_read) in reads.iter().enumerate() {
        if is_opp_orientation(&read, other_read) && is_overlap(&read, other_read) {
            let merge_result = attempt_merge(&read, other_read, min_overlap_bp);

            match merge_result {
                MergeResult::Discordant(_) => {
                    reads.remove(i);
                }
                MergeResult::NoMerge(_) => {}
                MergeResult::Merge(_) => {
                    reads.remove(i);
                }
            }
            return merge_result;
        }
    }

    MergeResult::NoMerge(())
}

pub fn construct_sequence<'a>(mut read_blueprint: IndexMap<i64, u8>) -> (i64, Vec<u8>) {
    let mut new_seq = Vec::new();

    read_blueprint.sort_unstable_keys();

    let start = read_blueprint
        .keys()
        .min()
        .expect("unable to find minimum genome pos");

    for base in read_blueprint.values() {
        new_seq.push(*base);
    }

    (*start, new_seq)
}

pub fn construct_read(
    original_read: &BamRecord,
    new_seq: Vec<u8>,
    mapper: &mut ReMapper,
    ref_seq: &Vec<u8>,
) -> BamRecord {
    let mut new_rec = original_read.clone();

    let (start, _end, cigar) = align_to_ref(mapper, &new_seq, &ref_seq);

    let qname = [new_rec.qname(), b":MERGED"].concat();
    new_rec.set(
        &qname,
        Some(&cigar),
        new_seq.as_slice(),
        vec![255; new_seq.len() as usize].as_slice(),
    );

    new_rec.set_pos(start as i64);
    return new_rec;
}

// with two overlapping reads, attempt to merge the reads
// halt if the reads have discordant sequence
pub fn attempt_merge(read_a: &BamRecord, read_b: &BamRecord, min_overlap_bp: usize) -> MergeResult {
    // check that these reads have opposing orientation
    let mut ra: IndexMap<i64, u8> = IndexMap::new();
    let mut rb: IndexMap<i64, u8> = IndexMap::new();

    let ras = read_a.seq().as_bytes();
    let rbs = read_b.seq().as_bytes();

    read_a.aligned_pairs().for_each(|pair| {
        ra.entry(pair[1]).or_insert(ras[pair[0] as usize]);
    });

    read_b.aligned_pairs().for_each(|pair| {
        rb.entry(pair[1]).or_insert(rbs[pair[0] as usize]);
    });

    let mut num_overlap = 0;
    let mut discordant = false;

    for (gpos, nuc) in ra {
        if let Some(other_nuc) = rb.get(&gpos) {
            if *other_nuc != nuc {
                warn!(
                    "\rDiscordant read detected! base a: {}, base b: {}",
                    str::from_utf8(&[nuc]).unwrap(),
                    str::from_utf8(&[*other_nuc]).unwrap(),
                );
                discordant = true;
                break;
            } else {
                num_overlap += 1;
            }
        } else {
            rb.entry(gpos).or_insert(nuc);
        }
    }
    if !discordant {
        if num_overlap >= min_overlap_bp {
            MergeResult::Merge(Some(rb))
        } else {
            MergeResult::NoMerge(())
        }
    } else {
        MergeResult::Discordant(())
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::record::record::BamRecord;
    use bio::alignment::pairwise::banded::Aligner;
    use bio::scores::blosum62;
    use crossbeam::channel::{bounded, Receiver, Sender};
    use rust_htslib::bam::{
        record::{Cigar, CigarString},
        Record,
    };

    fn create_bam_record(qname: &str, tid: i32, pos: i64, seq: &str, is_reverse: bool) -> Record {
        let mut record = Record::new();
        record.set_qname(qname.as_bytes());
        record.set_tid(tid);
        record.set_pos(pos);
        record.set(
            b"test",
            Some(&CigarString(vec![Cigar::Match(seq.len() as u32)])),
            seq.as_bytes(),
            vec![255; seq.len()].as_slice(),
        );
        if is_reverse {
            record.set_reverse();
            record.set_flags(0x2 | 0x10)
        } else {
            record.set_flags(0x2);
        }
        // record.set_flag(record.flag() & !0x4);
        println!("{}", record.is_reverse());
        record
    }

    #[test]
    fn test_is_opp_orientation() {
        let read_a = create_bam_record("read_a", 0, 10, "ATCG", false);
        let read_b = create_bam_record("read_b", 2, 10, "ATCG", true);

        assert!(is_opp_orientation(&read_a, &read_b));
        assert!(!is_opp_orientation(&read_a, &read_a));
    }

    #[test]
    fn test_is_overlap() {
        let read_a = create_bam_record("read_a", 0, 10, "ATCG", false);
        let read_b = create_bam_record("read_b", 0, 12, "ATCG", false);

        assert!(is_overlap(&read_a, &read_b));

        let read_c = create_bam_record("read_c", 0, 20, "ATCG", false);
        assert!(!is_overlap(&read_a, &read_c));
    }

    #[test]
    fn test_construct_sequence() {
        let mut read_blueprint = IndexMap::new();
        read_blueprint.insert(10, b'A');
        read_blueprint.insert(11, b'T');
        read_blueprint.insert(12, b'C');
        read_blueprint.insert(13, b'G');

        let (start, seq) = construct_sequence(read_blueprint);
        assert_eq!(start, 10);
        assert_eq!(seq, vec![b'A', b'T', b'C', b'G']);
    }

    // Test the handle_dupes function (simplified test)
    #[test]
    fn test_handle_dupes() {
        let mut umis_reads: IndexMap<String, ReadsAndCount<BamRecord>> = IndexMap::new();
        let mut merge_report = MergeReport::new();
        let (s, r): (Sender<Option<Record>>, Receiver<Option<Record>>) = bounded(10);
        umis_reads.insert(
            "umi1".to_string(),
            ReadsAndCount {
                reads: vec![
                    create_bam_record("read1", 0, 10, "ATCG", false),
                    create_bam_record("read2", 0, 13, "GATC", true),
                ],

                count: 2,
            },
        );

        let mapper: ReMapper = Aligner::new(-5, -1, blosum62, 19, 70);
        let ref_fasta = vec![b'A', b'T', b'C', b'G', b'A', b'T', b'C'];

        let results = handle_dupes(&mut umis_reads, mapper, &ref_fasta, 1, s.clone());
        for res in results {
            merge_report.count(res);
        }
        let out_read = r.recv().unwrap().unwrap();
        println!("{}", merge_report);
        println!("{:?}", out_read.cigar().to_string());
        println!("{:?}", out_read.seq().as_bytes());
        assert_eq!(out_read.seq().as_bytes(), b"ATCGATC");
        assert_eq!(out_read.cigar().to_string(), "7=");
    }
}
