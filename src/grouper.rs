use crate::group::bktree::NGramBKTree;
use crate::ngram::ngram;
use crate::GroupingMethod;
use indexmap::IndexSet;
use std::collections::{HashMap, VecDeque};
use std::sync::Arc;

pub struct Grouper<'a> {
    pub umis: &'a Vec<String>,
    pub ngram_maker: ngram::NgramMaker,
}

impl<'a> Grouper<'a> {
    pub fn new(umis: &'a Vec<String>, max_edit: usize, umi_len: usize) -> Self {
        let ngram_maker = ngram::NgramMaker::new(max_edit + 1, umi_len);

        Self {
            umis,
            ngram_maker,
        }
    }

    /// for every queried node, cluster its immediate offshoots distant by k edits, then retrieve the
    /// offshoots of each offshoot, until no more offshoots are found.
    pub fn visit_and_remove_all(
        &self,
        umi: &str,
        k: u8,
        counts: &HashMap<&'a str, i32>,
        bktree: &mut NGramBKTree,
    ) -> IndexSet<String> {
        let mut to_cluster = VecDeque::from([umi.to_string()]);
        let mut out = IndexSet::new();

        while let Some(root) = to_cluster.pop_front() {
            let max_count = (0.5 * (counts.get(root.as_str()).unwrap() + 1) as f32) as i32;
            let immediate = bktree.remove_near(root.as_str(), k, max_count, &self.ngram_maker);

            for c in immediate.iter().filter(|c| **c != root) {
                to_cluster.push_back(c.to_string());
            }

            out.extend(immediate);
        }

        out
    }

    pub fn cluster(
        &'a self,
        counts: HashMap<&'a str, i32>,
        grouping_method: Arc<&GroupingMethod>,
    ) -> impl Iterator<Item = IndexSet<String>> + use <'a> {
        let mut bktree = NGramBKTree::init_empty(None);
        bktree.count_map = counts.clone();

        for u in self.umis {
            bktree.populate_single(u, *counts.get(u.as_str()).unwrap(), &self.ngram_maker);
        }

        // cluster each umi based on the tree, only returning non-empty clusters.
        self.umis
            .iter()
            .map(move |u| self.visit_and_remove_all(u, 1, &counts, &mut bktree))
            .filter_map(|o| match o.is_empty() {
                true => None,
                false => Some(o),
            })
    }
}
