#![allow(dead_code)]

use crate::grouper::edit_distance;
use crate::ngram::ngram;
use indexmap::{IndexMap, IndexSet};
use std::cmp;
use std::collections::HashMap;
use std::ops::AddAssign;

/// Represents a Node of a BK-tree; each node contains its given string,
/// a collection of children strings (one per unique edit distance), and whether or not the node
/// has been invalidated. The minimum frequency of all child strings is also tracked.
#[derive(Debug)]
pub struct Node {
    umi: String,
    children: HashMap<usize, Node>,
    k: u32,
    exists: bool,
    count: usize,
    min_count: usize,
}

impl std::fmt::Display for Node {
    fn fmt(&self, _formatter: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        println!("Node: {}", self.umi);
        print!("Children:");
        for (k, node) in &self.children {
            println!("K: {k};\t{:?}", node)
        }

        println!("--------------------------");

        Ok(())
    }
}

impl Default for Node {
    fn default() -> Node {
        Node {
            umi: "".to_string(),
            children: HashMap::new(),
            k: 0,
            exists: true,
            count: 0,
            min_count: 0,
        }
    }
}

impl Node {
    pub fn new(umi: &str, count: usize) -> Self {
        let mut s = Self::default();
        s.umi = umi.to_string();
        s.count = count;
        s
    }

    /// attempt to add a string as a child to a node, unless a child with the same edit distance exists for
    /// the given node. Return the child if it exists.
    pub fn append_child(&mut self, k: usize, umi: &str, count: usize) -> Option<&mut Node> {
        // todo: rename
        if self.children.get(&k).is_none() {
            let c = Node::new(umi, count);

            self.children.entry(k).insert_entry(c);
            None
        } else {
            Some(self.children.get_mut(&k).unwrap())
        }
    }
}

/// Represents a collection of BK-trees, one for each unique ngram of an alphabet.
/// Each BK-tree is represented as its root node; see [Node] for more information.
pub struct NGramBKTree<'a> {
    ngram_tree_map: HashMap<String, Node>,
    count_map: HashMap<&'a str, usize>,
    root: Node,
    ngram_maker: ngram::NgramMaker,
    capacity: usize,
    ngrams: Vec<ngram::Ngram<'a>>,
}

impl std::fmt::Display for NGramBKTree<'_> {
    fn fmt(&self, _formatter: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        for (ngram, bktree) in &self.ngram_tree_map {
            println!("ngram: {ngram}\n{bktree}")
        }
        Ok(())
    }
}

impl<'a> NGramBKTree<'a> {
    pub fn init_empty(cap: Option<usize>, str_len: usize) -> Self {
        Self {
            ngram_tree_map: HashMap::with_capacity(cap.unwrap_or(100)),
            count_map: HashMap::with_capacity(cap.unwrap_or(100)),
            root: Node::default(),
            ngram_maker: ngram::NgramMaker::new(4, str_len),
            capacity: cap.unwrap_or(100),
            ngrams: vec!["NILL"; 4],
        }
    }

    pub fn is_empty(&self) -> bool {
        self.ngram_tree_map.is_empty()
    }

    /// Traverses a B-tree node-by-node until finding a node without a child at the given edit
    /// distance. Creates child node from query string and count once found.
    pub fn insert_raw_string(&mut self, node: &mut Node, umi: &str, count: usize) {
        let mut k: usize;
        let mut curr = node;

        k = edit_distance(&umi, &curr.umi);
        curr.count = cmp::min(curr.min_count, count);

        while let Some(node) = curr.append_child(k, umi, count) {
            curr = node;
            k = edit_distance(&umi, &curr.umi);
            curr.min_count = cmp::min(curr.min_count, count);
        }
    }

    pub fn remove_near(&mut self, umi: &str, k: usize, count: usize, max_edit: usize) {
        let mut found: IndexSet<ngram::Ngram> = IndexSet::new();
        let mut ngrams = std::mem::take(&mut self.ngrams);
        let mut ngram_map = std::mem::take(&mut self.ngram_tree_map);

        self.ngram_maker.ngrams_to_vec(umi, &mut ngrams);

        for n in ngrams {
            if let Some(curr) = ngram_map.get_mut(n) {
                self.recursive_remove_near(curr, umi, k, count, max_edit, &mut found);
            }
        }
    }

    pub fn recursive_remove_near(
        &mut self,
        curr: &'a mut Node,
        umi: &str,
        k: usize,
        max_count: usize,
        max_edit: usize,
        output_set: &mut IndexSet<&'a str>,
    ) {
        let dist = edit_distance(&curr.umi, umi);
        if dist <= k && curr.count <= max_count {
            output_set.insert(&curr.umi);
            curr.exists = false;
        }

        let mut min_count = self.count_map.get_mut(curr.umi.as_str()).unwrap();

        if !curr.children.is_empty() {
            let lo = cmp::max(dist - k, 0);
            let length = curr.children.len();
            let hi = cmp::min(dist + k, length - 1);

            for i in 0..length {
                if let Some(child) = curr.children.get_mut(&i) {
                    if i >= lo && i <= hi && child.min_count <= max_count {
                        self.recursive_remove_near(child, umi, k, max_count, max_edit, output_set);
                    }

                    min_count = cmp::min(min_count, &mut child.min_count);
                }
            }
        }
        unimplemented!();
    }

    /// given a hashmap relating containing all strings related to a given ngram, create a BK-tree
    /// for each ngram holding all associated strings.
    pub fn populate_from_ngram_map(&mut self, ngram_map: IndexMap<ngram::Ngram, IndexSet<&str>>) {
        let mut ngram_tree_map = std::mem::take(&mut self.ngram_tree_map);
        let mut k;

        for (ngram, mut neighbors) in ngram_map {
            let en = ngram_tree_map.get(ngram);
            let mut neighbors = neighbors.drain(..);

            if en.is_none() {
                let first_umi = neighbors.next().unwrap();

                ngram_tree_map.insert(
                    ngram.to_string(),
                    Node::new(&first_umi, *self.count_map.get(first_umi).unwrap()),
                );
            }

            let en = ngram_tree_map.get_mut(ngram).unwrap();

            for n in neighbors {
                k = edit_distance(&en.umi, &n);
                en.append_child(k, &n, *self.count_map.get(&n).unwrap());
            }
        }

        self.ngram_tree_map = ngram_tree_map;
    }
}

#[cfg(test)]
use std::iter;
#[test]
fn init_test() {
    let cap = 42069;
    let ngram_bk = NGramBKTree::init_empty(Some(cap), 12);
    assert!(ngram_bk.count_map.capacity() >= cap);
    assert!(ngram_bk.ngram_tree_map.capacity() >= cap);
}

#[test]
fn populate_from_umis() {
    let umi_a = "GATACAGA";
    let umi_b = "GATACAGT";
    let umi_c = "TATACATA";

    let umis = iter::repeat(umi_a)
        .take(5)
        .chain(iter::repeat(umi_b).take(3))
        .chain(iter::repeat(umi_c).take(1))
        .collect::<Vec<&str>>();

    let mut counts: HashMap<&str, usize> = HashMap::new();
    let mut ngram_map: IndexMap<ngram::Ngram, IndexSet<&str>> = IndexMap::new();
    let mut ngram_vec = vec!["NULL"; 2];
    let mut ngram_maker = ngram::NgramMaker::new(2, umi_a.len());

    for u in umis {
        counts.entry(u).or_insert(0).add_assign(1);
        ngram_maker.ngrams_to_vec(u, &mut ngram_vec);
        for ngram in &ngram_vec {
            ngram_map.entry(ngram).or_insert(IndexSet::new()).insert(u);
        }
    }

    let mut bktree = NGramBKTree::init_empty(None, umi_a.len());
    bktree.count_map = counts;
    bktree.populate_from_ngram_map(ngram_map);
    println! {"{bktree}"};
}
