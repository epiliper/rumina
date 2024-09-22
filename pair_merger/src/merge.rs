use indexmap::IndexMap;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use std::collections::HashMap;

pub fn handle_dupes(umis_reads: HashMap<String, Vec<Record>>) -> Vec<Record> {
    let mut corrected_reads: Vec<Record> = Vec::new();

    for (_umi, reads) in umis_reads {
        match reads.len() {
            1 => panic!("Error: 1 read found for UMI marked as duplicate. Has the file been modified? Exiting..."),
            2 => {
                let read_a = reads.first().unwrap();
                let read_b = reads.last().unwrap();

                if let Some(merged_seq) = get_overlap(read_a, read_b) {
                    let new_seq = construct_sequence(merged_seq);
                    let new_record = construct_read(read_a, new_seq);
                    corrected_reads.push(new_record);
                }
                todo!();
            }
            _ => println!("Warning: 3 or more duplicate UMIs detected. Deduplicating the first two"),
        }
    }
    todo!();
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
    new_rec.set(
        original_read.qname(),
        None,
        new_seq.as_slice(),
        vec![255, new_seq.len() as u8].as_slice(),
    );
    return new_rec;
}

pub fn get_overlap(read_a: &Record, read_b: &Record) -> Option<IndexMap<i64, u8>> {
    // check that these reads have opposing orientation
    if (read_a.is_reverse() && !read_b.is_reverse()) | (!read_a.is_reverse() && read_b.is_reverse())
    {
        print!("Forward and reverse read pair identified.");

        let ra = read_a.aligned_pairs();
        let rb = read_a.aligned_pairs();

        // this will hold the reconstructed read
        // key: index of base along read
        // value: base
        let mut new_seq: IndexMap<i64, u8> = IndexMap::new();

        let mut discordant = false;

        let ras = read_a.seq().as_bytes();
        let rab = read_b.seq().as_bytes();

        for (a_seqi, b_seqi) in ra.zip(rb) {
            let base_a = ras[a_seqi[1] as usize];
            let base_b = rab[b_seqi[1] as usize];

            let genome_a = a_seqi[0];
            let genome_b = b_seqi[0];

            match genome_a == genome_b {
                true => {
                    // overlapping pos on genome

                    if base_a == base_b {
                        new_seq.insert(a_seqi[1], base_b); // overlap bases match
                    } else {
                        // if bases don't match at the same genome position, read pair is
                        // discordant.
                        println!("Discordant read pair found. Discarding...");
                        discordant = true;
                        break;
                    }
                }

                false => {
                    new_seq.insert(a_seqi[1], base_a);
                    new_seq.insert(b_seqi[1], base_b);
                }
            }
        }

        if discordant {
            return None;
        } else {
            return Some(new_seq);
        }
    }
    return None;
}
