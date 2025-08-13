use smol_str::SmolStr;
use std::hash::{Hash, Hasher};

use crate::record::Record;
use indexmap::IndexMap;

// every position has a key
pub type PositionKey<T> = IndexMap<i64, KeyUMI<T>>;

// every key has a set of UMIs and associated reads
pub type KeyUMI<T> = IndexMap<u64, UmiReadMap<T>>;

#[derive(Debug)]
pub struct ReadsAndCount<T>
where
    T: Record,
{
    pub reads: Vec<T>,
    pub count: i32,
}

// when a read is added to ReadsAndCount, increase the read count at x umi
impl<T: Record> ReadsAndCount<T> {
    pub fn up(&mut self, read: T) {
        self.count += 1;
        self.reads.push(read);
    }
}

/// Maps UMIs to associated reads, where reads are stratified further by sequence (see [SeqEntry]).
pub type UmiReadMap<T> = IndexMap<SmolStr, (i32, SeqMap<T>)>;

/// Associates all reads sharing a given sequence.
pub type SeqMap<T> = IndexMap<u64, SeqEntry<T>>;

pub trait ReadStore<T: Record> {
    fn combine(&mut self, other: SeqMap<T>, retain_all: bool);
    fn intake(&mut self, read: T, retain_all: bool) -> u8;
}

impl<T: Record> ReadStore<T> for SeqMap<T> {
    /// Combine two [SeqMap]s into one
    fn combine(&mut self, mut other: SeqMap<T>, retain_all: bool) {
        other.drain(..).for_each(|(other_seq, mut seq_entry)| {
            let mut e = self.entry(other_seq).or_insert(SeqEntry::new(retain_all));
            seq_entry
                .reads
                .drain(..)
                .for_each(|read| {((e.up_method)(&mut e, read));})
            // if let Some(mut seq) = self.get_mut(&other_seq) {
            //     seq_entry.reads.drain(..).for_each(|read| {
            //         (seq.up_method)(&mut seq, read);
            //     })
        });
    }

    /// Update [Self] with a new read
    fn intake(&mut self, read: T, retain_all: bool) -> u8 {
        let mut h = std::hash::DefaultHasher::new();
        read.seq_str().hash(&mut h);
        let mut e = self.entry(h.finish()).or_insert(SeqEntry::new(retain_all));

        e.count += 1;
        (e.up_method)(&mut e, read)
    }
}

#[derive(Debug)]
/// Holds a list of records, intended to store all reads matching in sequence (see [SeqMap]).
/// [Self::count] tracks all reads with matching sequence, regardless of whether all or one is retained.
///
/// Can hold all reads per sequence, or one read per sequence, in which case the read with the best
/// sequence quality is retained per sequence.
pub struct SeqEntry<T>
where
    T: Record,
{
    // making this an array for iteration
    pub reads: Vec<T>,
    pub count: i32,
    pub qual_sum: u32,
    pub up_method: for<'a> fn(&'a mut Self, read: T) -> u8,
}

impl<T: Record> SeqEntry<T> {
    /// only keep one read per sequence
    pub fn up_keep_single(&mut self, read: T) -> u8 {
        let s: u32 = read.qual().iter().map(|a| u32::try_from(*a).unwrap()).sum();

        let ret = match self.reads.is_empty() {
            true => 1,
            false => 0,
        };

        // using >= here in case a read is Q of 0 (though such crap reads should be filtered out
        // first)
        if s >= self.qual_sum {
            self.reads = vec![read];
            self.qual_sum = s;
        }

        assert_eq!(self.reads.len(), 1);
        self.count += 1;

        ret
    }

    /// keep all
    pub fn up_group(&mut self, read: T) -> u8 {
        self.reads.push(read);
        self.count += 1;
        1
    }

    pub fn new(group: bool) -> Self {
        // determine whether to retain a single read per sequence, or all reads
        let up_method = match group {
            true => Self::up_group,
            false => Self::up_keep_single,
        };

        let reads = vec![];
        let count = 0;
        let qual_sum = 0;

        Self {
            reads,
            count,
            qual_sum,
            up_method,
        }
    }
}
