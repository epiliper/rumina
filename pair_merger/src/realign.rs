use bio::alignment::pairwise::banded::*;
use bio::io::fasta;
use bio::scores::blosum62;
use rust_htslib::bam::record::{Cigar, CigarString};

pub type ReMapper = Aligner<fn(u8, u8) -> i32>;

pub fn init_remapper(ref_fasta_file: String) -> (ReMapper, Vec<u8>) {
    let aligner: ReMapper = Aligner::new(-5, -1, blosum62, 19, 70);

    let reader = fasta::Reader::from_file(ref_fasta_file).expect("Error reading fasta file!");
    let reference = reader
        .records()
        .next()
        .unwrap()
        .expect("Error retrieving reference genome from fasta file!");

    let ref_seq = reference.seq().to_owned();
    (aligner, ref_seq)
}

pub fn align_to_ref(
    aligner: &mut ReMapper,
    record_seq: &Vec<u8>,
    ref_seq: &Vec<u8>,
) -> (usize, usize, CigarString) {
    let aln = aligner.semiglobal(record_seq, ref_seq);

    (aln.ystart, aln.yend, parse_cigar_string(&aln.cigar(false)))
}

fn parse_cigar_string(cigar_str: &str) -> CigarString {
    let mut ops = Vec::new();
    let mut current_num = 0;
    for c in cigar_str.chars() {
        match c {
            '0'..='9' => current_num = current_num * 10 + (c as u32 - '0' as u32),
            _ => {
                let op = match c {
                    'M' => Cigar::Match(current_num),
                    'I' => Cigar::Ins(current_num),
                    'D' => Cigar::Del(current_num),
                    'N' => Cigar::RefSkip(current_num),
                    'S' => Cigar::SoftClip(current_num),
                    'H' => Cigar::HardClip(current_num),
                    'P' => Cigar::Pad(current_num),
                    'X' => Cigar::Diff(current_num),
                    '=' => Cigar::Equal(current_num),
                    _ => panic!("Invalid CIGAR operation: {}", c),
                };
                ops.push(op);
                current_num = 0;
            }
        }
    }

    CigarString(ops)
}
