use crate::record::Record;
use indexmap::IndexMap;

pub type PositionKey<T> = IndexMap<i64, KeyUMI<T>>; //every position has a key
pub type KeyUMI<T> = IndexMap<u64, UmiReadMap<T>>;
// pub type UMIReads<T> = IndexMap<String, SeqEntry<T>>; // every UMI has a set of reads

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
pub type UmiReadMap<T> = IndexMap<String, (i32, SeqMap<T>)>;

// pub struct SeqMap<T: Record>(IndexMap<String, SeqEntry<T>>);
pub type SeqMap<T> = IndexMap<String, SeqEntry<T>>;

pub trait ReadStore<T: Record> {
    fn combine(&mut self, other: SeqMap<T>);
    fn intake(&mut self, read: T, retain_all: bool);
}

impl<T: Record> ReadStore<T> for SeqMap<T> {
    fn combine(&mut self, mut other: SeqMap<T>) {
        other.drain(..).for_each(|(other_seq, mut seq_entry)| {
            if let Some(mut seq) = self.get_mut(&other_seq) {
                seq_entry
                    .reads
                    .drain(..)
                    .for_each(|read| (seq.up_method)(&mut seq, read))
            };
        })
    }

    fn intake(&mut self, read: T, retain_all: bool) {
        let mut e = self.entry(read.seq()).or_insert(SeqEntry::new(retain_all));

        e.count += 1;
        (e.up_method)(&mut e, read);
    }
}

#[derive(Debug)]
/// Holds a list of records, intended to store all reads matching in sequence (see [UmiEntry]).
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
    pub up_method: for<'a> fn(&'a mut Self, read: T),
}

impl<T: Record> SeqEntry<T> {
    /// only keep one read per sequence
    pub fn up_keep_single(&mut self, read: T) {
        let s: u32 = read.qual().iter().map(|a| u32::try_from(*a).unwrap()).sum();
        if s > self.qual_sum {
            self.reads = vec![read];
            self.qual_sum = s;
        }

        assert_eq!(self.reads.len(), 1);
        self.count += 1;
    }

    /// keep all
    pub fn up_group(&mut self, read: T) {
        self.reads.push(read);
        self.count += 1;
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
