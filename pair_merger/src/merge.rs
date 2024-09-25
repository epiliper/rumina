use indexmap::IndexMap;
use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Cigar, CigarString},
    Record,
};
use std::collections::HashMap;

pub fn handle_dupes(umis_reads: &mut HashMap<String, Vec<Record>>) -> Vec<Record> {
    let mut corrected_reads: Vec<Record> = Vec::new();
    for (_umi, reads) in umis_reads {
        match reads.len() {
            1 => {
                println!("Warning: 1 read found for UMI marked as duplicate. Has the file been modified?");
                // todo, remove this read, don't add.
                corrected_reads.extend(reads.drain(..));
            }
            2 => {
                let read_a = reads.first().unwrap();
                let read_b = reads.last().unwrap();

                if let Some(merged_seq) = get_overlap(read_a, read_b) {
                    let (start_pos, new_seq) = construct_sequence(merged_seq);
                    let new_record = construct_read(read_a, start_pos, new_seq);
                    corrected_reads.push(new_record);
                }
            }
            _ => {
                println!("Warning: 3 or more duplicate UMIs detected. Deduplicating the first two")
            }
        }
    }

    corrected_reads
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

pub fn construct_read(original_read: &Record, start_pos: i64, new_seq: Vec<u8>) -> Record {
    let mut new_rec = original_read.clone();
    let qname = [new_rec.qname(), b":MERGED"].concat();
    new_rec.set(
        &qname,
        // None,
        Some(&CigarString(vec![Cigar::Match(new_seq.len() as u32)])),
        new_seq.as_slice(),
        vec![255; new_seq.len() as usize].as_slice(),
    );

    new_rec.set_pos(start_pos);
    return new_rec;
}

pub fn get_overlap(read_a: &Record, read_b: &Record) -> Option<IndexMap<i64, u8>> {
    // check that these reads have opposing orientation
    if (read_a.is_reverse() && !read_b.is_reverse()) | (!read_a.is_reverse() && read_b.is_reverse())
    {
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

        let mut discordant = false;
        let mut overlap = false;

        for (gpos, nuc) in ra {
            if let Some(other_nuc) = rb.get(&gpos) {
                if *other_nuc != nuc {
                    println!("Discordant read detected! base a: {nuc}, base b: {other_nuc}");
                    discordant = true;
                    break;
                } else {
                    overlap = true;
                }
            } else {
                rb.entry(gpos).or_insert(nuc);
            }
        }

        if !discordant && overlap {
            return Some(rb);
        } else {
            return None;
        }
    } else {
        None
    }
}
