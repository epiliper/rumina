// #![allow(dead_code)]

use crate::ngram::ngram;
use indexmap::IndexSet;
use std::cell::RefCell;
use std::cmp;
use std::collections::HashMap;
use std::collections::VecDeque;
use std::rc::Rc;

pub fn hamming(ua: &str, ub: &str) -> usize {
    strsim::hamming(ua, ub).unwrap()
}

/// Represents a Node of a BK-tree; each node contains its given string,
/// a collection of children strings (one per unique edit distance), and whether or not the node
/// has been invalidated. The minimum frequency of all child strings is also tracked.
#[derive(Debug)]
pub struct Node {
    umi: String,
    children: HashMap<u8, Rc<RefCell<Node>>>,
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
        println!("Node with count {}: {}", self.count, self.umi);
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
        s.min_count = count;
        s
    }

    /// attempt to add a string as a child to a node, unless a child with the same edit distance exists for
    /// the given node. Return the child if it exists.
    pub fn append_child(&mut self, k: u8, umi: &str, count: i32) -> Option<Rc<RefCell<Node>>> {
        if self.children.get(&k).is_none() {
            let c = Rc::new(RefCell::new(Node::new(umi, count)));

            self.children.entry(k).insert_entry(c);

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
    pub ngram_tree_map: HashMap<String, Rc<RefCell<Node>>>,
    pub count_map: HashMap<&'a str, i32>,
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
    pub fn init_empty(cap: Option<usize>) -> Self {
        Self {
            ngram_tree_map: HashMap::with_capacity(cap.unwrap_or(100)),
            count_map: HashMap::with_capacity(cap.unwrap_or(100)),
        }
    }

    pub fn contains(&self, s: &str) -> bool {
        self.count_map.contains_key(s)
    }

    /// Traverses a B-tree node-by-node until finding a node without a child at the given edit
    /// distance. Creates child node from query string and count once found.
    pub fn insert_raw_string(&self, node: Rc<RefCell<Node>>, umi: &str, count: i32) {
        let mut k;
        let mut child_found;
        let mut curr = node.clone();

        loop {
            {
                let mut curr = curr.borrow_mut();
                k = hamming(&umi, &curr.umi) as u8;
                if k == 0 {
                    // umi already in tree
                    break;
                }
                curr.min_count = cmp::min(curr.min_count, count);
                child_found = curr.append_child(k, umi, count);
            }

            if let Some(child) = child_found {
                curr = child;
            } else {
                break;
            }
        }
    }

    pub fn remove_near<'b>(
        &mut self,
        umi: &'b str,
        k: u8,
        max_count: i32,
        ngm: &ngram::NgramMaker,
    ) -> IndexSet<String> {
        let mut found = IndexSet::new();
        for ngram in ngm.ngrams(umi).iter() {
            if let Some(node) = self.ngram_tree_map.get(ngram) {
                self.remove_near_stack(node.clone(), umi, 0, i32::MAX, &mut found);
                self.remove_near_stack(node.clone(), umi, k, max_count, &mut found);
            }
        }

        self.prune(&found);
        found
    }

    pub fn prune(&mut self, to_remove: &IndexSet<String>) {
        for n in to_remove {
            self.count_map.remove(n.as_str());
        }
    }

    pub fn remove_near_stack(
        &self,
        node: Rc<RefCell<Node>>,
        umi: &str,
        k: u8,
        max_count: i32,
        output: &mut IndexSet<String>,
    ) {
        let mut visited: VecDeque<Rc<RefCell<Node>>> = VecDeque::from([node.clone()]);

        while let Some(node_ref) = visited.pop_front() {
            let node = node_ref.borrow();

            let dist = hamming(&node.umi, umi) as u8;
            let min_dist = dist.saturating_sub(k);
            let max_dist = dist.saturating_add(k);

            // check node children for within k threshold
            for i in min_dist..=max_dist {
                if let Some(child) = node.children.get(&i) {
                    visited.push_back(child.clone());
                }
            }

            // also add the current node if it's within k edits
            if dist <= k && node.count <= max_count && self.contains(node.umi.as_str()) {
                output.insert(node.umi.clone());
            }
        }
    }

    pub fn populate_single(&mut self, s: &str, count: i32, ngm: &ngram::NgramMaker) {
        for n in ngm.ngrams(s).iter() {
            if let Some(node) = self.ngram_tree_map.get(n) {
                self.insert_raw_string(node.clone(), s, count);
            } else {
                self.ngram_tree_map
                    .insert(n.to_string(), Rc::new(RefCell::new(Node::new(s, count))));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ngram::ngram;
    use std::collections::HashMap;

    #[test]
    fn test_ngram_edge_case_1() {
        // UMIs "ACG" and "GCA" are 1 edit distance apart, but share no 2-grams.
        let umi_a = "ACT";
        let umi_b = "ATT";

        let mut counts: HashMap<&str, i32> = HashMap::new();
        counts.insert(umi_a, 10);
        counts.insert(umi_b, 5);

        let ngram_maker = ngram::NgramMaker::new(2, umi_a.len());

        let mut bktree = NGramBKTree::init_empty(None);
        bktree.count_map = counts.clone();

        // Populate BK-tree directly to avoid IndexMap/neighbor issues
        for (umi, count) in &counts {
            for ngram_str in ngram_maker.ngrams(umi).iter() {
                let ngram = ngram_str.to_string();
                if !bktree.ngram_tree_map.contains_key(&ngram) {
                    bktree
                        .ngram_tree_map
                        .insert(ngram.clone(), Rc::new(RefCell::new(Node::new(umi, *count))));
                }
                let node = bktree.ngram_tree_map.get(&ngram).unwrap().clone();
                if *umi != node.borrow().umi {
                    bktree.insert_raw_string(node.clone(), umi, *count);
                }
            }
        }

        let res = bktree.remove_near(umi_a, 1, 10, &ngram_maker);

        // umi_b should be found, even though it shares no n-grams with umi_a's n-gram root.
        assert!(res.contains(umi_a), "Expected to find umi_a");
        assert!(res.contains(umi_b), "Expected to find umi_b");
    }

    #[test]
    fn test_ngram_edge_case_2() {
        // UMIs differing by a single insertion/deletion that disrupts all n-grams.
        let umi_a = "AAAA";
        let umi_b = "AAAT";
        let umi_c = "ATAT";

        let mut counts: HashMap<&str, i32> = HashMap::new();
        counts.insert(umi_a, 10);
        counts.insert(umi_b, 5);
        counts.insert(umi_c, 7);

        let ngram_maker = ngram::NgramMaker::new(2, umi_a.len());

        let mut bktree = NGramBKTree::init_empty(None);
        bktree.count_map = counts.clone();

        // Populate BK-tree directly (same as above)
        for (umi, count) in &counts {
            for ngram_str in ngram_maker.ngrams(umi).iter() {
                let ngram = ngram_str.to_string();
                if !bktree.ngram_tree_map.contains_key(&ngram) {
                    bktree
                        .ngram_tree_map
                        .insert(ngram.clone(), Rc::new(RefCell::new(Node::new(umi, *count))));
                }
                let node = bktree.ngram_tree_map.get(&ngram).unwrap().clone();
                if *umi != node.borrow().umi {
                    bktree.insert_raw_string(node.clone(), umi, *count);
                }
            }
        }

        let res = bktree.remove_near(umi_a, 1, 10, &ngram_maker);

        assert!(res.contains(umi_a), "Expected to find umi_a");
        assert!(res.contains(umi_b), "Expected to find umi_b");
        assert!(!res.contains(umi_c), "Did not expect to find umi_c");
    }
}
