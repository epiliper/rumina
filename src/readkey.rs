use std::hash::{DefaultHasher, Hash, Hasher};

// this module is responsible for creating a key for batching reads
// basically, readkey stores batching-relevant features of the read,
// which are then combined and hashed to make a key.

pub struct ReadKey {
    pub length: usize,
    pub reverse: bool,
    pub chr: usize,
}

impl Hash for ReadKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.length.hash(state);
        self.reverse.hash(state);
        self.chr.hash(state);
    }
}

impl PartialEq for ReadKey {
    fn eq(&self, other: &Self) -> bool {
        self.length == other.length && self.reverse == other.reverse
    }
}

impl ReadKey {
    pub fn get_key(self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }
}
