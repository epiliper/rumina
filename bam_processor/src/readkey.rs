use std::hash::{DefaultHasher, Hash, Hasher};

pub struct ReadKey {
    pub length: i64,
    pub reverse: bool,
}

impl Hash for ReadKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (self.length * self.reverse as i64).hash(state)
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
