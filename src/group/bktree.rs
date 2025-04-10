#![allow(dead_code)]

use crate::grouper::levenshtein;
use crate::ngram::ngram;
use indexmap::{IndexMap, IndexSet};
use std::cell::RefCell;
use std::cmp;
use std::collections::HashMap;
use std::collections::VecDeque;
use std::ops::AddAssign;
use std::rc::Rc;

/// Represents a Node of a BK-tree; each node contains its given string,
/// a collection of children strings (one per unique edit distance), and whether or not the node
/// has been invalidated. The minimum frequency of all child strings is also tracked.
#[derive(Debug)]
pub struct Node {
    umi: String,
    children: HashMap<i8, Rc<RefCell<Node>>>,
    k: i8,
    subtree_exists: bool,
    count: i32,
    min_count: i32,
}

impl std::hash::Hash for Node {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.umi.hash(state);
    }
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.umi == other.umi
    }
}

impl Eq for Node {}

impl std::fmt::Display for Node {
    fn fmt(&self, _formatter: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        println!("Node: {}", self.umi);
        println!("Children:");
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
            subtree_exists: false,
            count: 0,
            min_count: 0,
        }
    }
}

impl Node {
    pub fn new(umi: &str, count: i32) -> Self {
        let mut s = Self::default();
        s.umi = umi.to_string();
        s.count = count;
        s
    }

    pub fn has_subtree(&self, k: i8) -> bool {
        self.children.contains_key(&k)
    }

    /// attempt to add a string as a child to a node, unless a child with the same edit distance exists for
    /// the given node. Return the child if it exists.
    pub fn append_child(&mut self, k: i8, umi: &str, count: i32) -> Option<Rc<RefCell<Node>>> {
        if self.children.get(&k).is_none() {
            let c = Rc::new(RefCell::new(Node::new(umi, count)));

            self.children.entry(k).insert_entry(c);

            self.children
                .entry(k)
                .and_modify(|k| k.borrow_mut().subtree_exists = true);

            self.subtree_exists = true;

            None
        } else {
            Some(self.children.get(&k).unwrap().clone())
        }
    }
}

/// Represents a collection of BK-trees, one for each unique ngram of an alphabet.
/// Each BK-tree is represented as its root node; see [Node] for more information.
pub struct NGramBKTree<'a> {
    ngram_tree_map: HashMap<String, Rc<RefCell<Node>>>,
    pub count_map: HashMap<&'a str, i32>,
    root: Node,
    ngram_maker: ngram::NgramMaker<'a>,
    capacity: usize,
    ngrams: Vec<ngram::Ngram<'a>>,
    str_len: usize,
}

impl std::fmt::Display for NGramBKTree<'_> {
    fn fmt(&self, _formatter: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        for (ngram, bktree) in &self.ngram_tree_map {
            println!("ngram: {ngram}\n{}", bktree.borrow())
        }
        Ok(())
    }
}

impl<'a> NGramBKTree<'a> {
    pub fn init_empty(cap: Option<usize>, max_edit: usize, str_len: usize) -> Self {
        Self {
            ngram_tree_map: HashMap::with_capacity(cap.unwrap_or(100)),
            count_map: HashMap::with_capacity(cap.unwrap_or(100)),
            root: Node::default(),
            ngram_maker: ngram::NgramMaker::new(max_edit + 1, str_len),
            capacity: cap.unwrap_or(100),
            ngrams: vec!["NILL"; 4],
            str_len,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.ngram_tree_map.is_empty()
    }

    /// Traverses a B-tree node-by-node until finding a node without a child at the given edit
    /// distance. Creates child node from query string and count once found.
    pub fn insert_raw_string(&mut self, node: Rc<RefCell<Node>>, umi: &str, count: i32) {
        let mut k;
        let mut curr = node.borrow_mut();

        k = levenshtein(&umi, &curr.umi) as i8;
        curr.count = cmp::min(curr.min_count, count);

        while let Some(child) = curr.append_child(k, umi, count) {
            k = levenshtein(&child.borrow().umi, umi) as i8;
            // update count
        }
    }

    pub fn remove_near(&mut self, umi: &'a str, k: i8, max_count: i32) -> IndexSet<String> {
        let mut ngrams = std::mem::take(&mut self.ngrams);
        let mut found: IndexSet<String> = IndexSet::new();
        self.ngram_maker.ngrams_to_vec(umi, &mut ngrams);

        for ngram in ngrams {
            if let Some(node) = self.ngram_tree_map.get(ngram) {
                self.remove_near_stack(node.clone(), umi, k, max_count, &mut found);
            }
        }

        return found;
    }

    pub fn remove_near_stack(
        &self,
        node: Rc<RefCell<Node>>,
        umi: &str,
        k: i8,
        max_count: i32,
        output: &mut IndexSet<String>,
    ) {
        let mut visited: VecDeque<Rc<RefCell<Node>>> = VecDeque::from([node.clone()]);

        while !visited.is_empty() {
            let node_ref = visited.pop_front().unwrap();
            let node = node_ref.borrow_mut();

            let dist = levenshtein(&node.umi, umi) as i8;
            let exists = self.count_map.contains_key(node.umi.as_str());

            if dist <= k && exists && node.count <= max_count {
                output.insert(node.umi.clone());
            }

            if node.subtree_exists {
                let low = cmp::max(dist - k, 0);
                let length = self.str_len as i8;
                let high = cmp::min(dist + k, length - 1);

                for i in low..=high {
                    if let Some(child) = node.children.get(&i) {
                        visited.push_back(child.clone());
                    }
                }
            }
        }
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
                    Rc::new(RefCell::new(Node::new(
                        &first_umi,
                        *self.count_map.get(first_umi).unwrap(),
                    ))),
                );
            }

            let mut en = ngram_tree_map.get_mut(ngram).unwrap().borrow_mut();

            for n in neighbors {
                k = levenshtein(&en.umi, &n) as i8;
                en.append_child(k, &n, *self.count_map.get(&n).unwrap() as i32);
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
    let ngram_bk = NGramBKTree::init_empty(Some(cap), 12, 1);
    assert!(ngram_bk.count_map.capacity() >= cap);
    assert!(ngram_bk.ngram_tree_map.capacity() >= cap);
}

#[test]
fn populate_from_umis() {
    let umi_a = "GATACAGA";
    let umi_b = "GTTACAGT";
    let umi_c = "GATACAGT";
    let umi_d = "TATACATA";

    let umis = iter::repeat(umi_a)
        .take(5)
        .chain(iter::repeat(umi_b).take(3))
        .chain(iter::repeat(umi_c).take(3))
        .chain(iter::repeat(umi_d).take(1))
        .collect::<Vec<&str>>();

    let mut counts: HashMap<&str, i32> = HashMap::new();
    let mut ngram_map: IndexMap<ngram::Ngram, IndexSet<&str>> = IndexMap::new();
    let mut ngram_maker = ngram::NgramMaker::new(2, umi_a.len());
    let mut ngram_vec = vec!["NULL"; ngram_maker.num_chunks];

    for u in umis {
        counts.entry(u).or_insert(0).add_assign(1);
        ngram_maker.ngrams_to_vec(u, &mut ngram_vec);
        for ngram in &ngram_vec {
            ngram_map.entry(ngram).or_insert(IndexSet::new()).insert(u);
        }
    }

    let mut bktree = NGramBKTree::init_empty(None, 1, umi_a.len());
    bktree.count_map = counts;
    bktree.ngrams = ngram_vec.clone();
    bktree.ngram_maker = ngram_maker;

    bktree.populate_from_ngram_map(ngram_map);
    let res = bktree.remove_near(umi_a, 1, 5);
    println! {"found for {umi_a}: {:?}", res};

    bktree.ngrams = ngram_vec;

    let res = bktree.remove_near(umi_b, 1, 5);
    println! {"found for {umi_b}: {:?}", res};
}
