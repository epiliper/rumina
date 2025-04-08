use indexmap::IndexMap;
use std::collections::HashMap;

pub struct Node {
    children: Vec<String>,
    k: u32,
    exists: bool,
    count: u64,
}

pub struct NGramBKTree {
    ngram_tree_map: HashMap<String, Node>,
    count_map: HashMap<String, i64>,
}

impl NGramBKTree {
    pub fn init_empty(cap: Option<usize>) -> Self {
        Self {
            ngram_tree_map: HashMap::with_capacity(cap.unwrap_or(100)),
            count_map: HashMap::with_capacity(cap.unwrap_or(100)),
        }
    }

    pub fn insert(&mut self, umi: String, count: u64) {}

    pub fn load_ngrams(&mut self) {
        unimplemented!();
    }
}
