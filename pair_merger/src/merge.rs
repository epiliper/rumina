use indexmap::IndexMap;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use std::collections::HashMap;
use std::str;

pub fn handle_dupes(umis_reads: &HashMap<String, Vec<Record>>) -> Vec<Record> {
    let mut corrected_reads: Vec<Record> = Vec::new();
    for (_umi, reads) in umis_reads {
        match reads.len() {
            1 => panic!("Error: 1 read found for UMI marked as duplicate. Has the file been modified? Exiting..."),
            2 => {
                let read_a = reads.first().unwrap();
                let read_b = reads.last().unwrap();

                // println!("{}\n{}",
                //     str::from_utf8(&read_a.seq().as_bytes()).unwrap(),
                //     str::from_utf8(&read_b.seq().as_bytes()).unwrap(),
                //     );

                if let Some(merged_seq) = get_overlap(read_a, read_b) {
                    let new_seq = construct_sequence(merged_seq);
                    let new_record = construct_read(read_a, new_seq);
                    // println!("{}\n----------------------", 
                    //     str::from_utf8(&new_record.seq().as_bytes()).unwrap(),
                    //     );
                    corrected_reads.push(new_record);
                    // corrected_reads.push(read_a);
                    // corrected_reads.push(read_b);
                }

                else {
                    println!("-----------------------")
                }
            }
            _ => println!("Warning: 3 or more duplicate UMIs detected. Deduplicating the first two"),
        }
    }

    corrected_reads
}

pub fn construct_sequence<'a>(mut read_blueprint: IndexMap<i64, u8>) -> Vec<u8> {
    let mut new_seq = Vec::new();
    read_blueprint.sort_unstable_keys();
    for base in read_blueprint.values() {
        new_seq.push(*base);
    }
    new_seq
}

pub fn construct_read(original_read: &Record, new_seq: Vec<u8>) -> Record {
    let mut new_rec = original_read.clone();
    let qname = [new_rec.qname(), b":MERGED"].concat();
    new_rec.set(
        &qname,
        // Some(&original_read.cigar()),
        None,
        new_seq.as_slice(),
        vec![255; new_seq.len() as usize].as_slice(),
    );
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
