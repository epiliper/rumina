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

    /// Represents the main functionality of [Grouper]: for all UMIs loaded into the struct,
    /// split each into ngrams and load them into BK trees according to the user-specified grouping
    /// method.
    ///
    /// Once all UMIs have been loaded, cluster each one with their predicted offshoots, which are
    /// determined by these parameters:
    ///
    /// k: the maximum edit distance between two UMIs for them to be clustered together
    ///
    /// percentage: the maximum percentage of a UMI's count another UMI must have to be considered
    /// an offshoot.
    ///
    /// For example: 
    /// k = 1, percentage = 0.5
    /// UMI a: ATCG, count: 10
    /// UMI b: ATCC, count: 4
    /// UMI c: TTCC, count: 1
    /// UMI d: ATGA, count: 40
    ///
    /// { a -> b -> c } { d }
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

    /// for every queried node, cluster its immediate offshoots distant by <= k edits and the given count percentage, then retrieve the
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
            // TODO: make percentage adjustable
            let max_count = (0.5 * (counts.get(root.as_str()).unwrap() + 1) as f32) as i32;
            let immediate = bktree.remove_near(root.as_str(), k, max_count, &self.ngram_maker);

            for c in immediate.iter().filter(|c| **c != root) {
                to_cluster.push_back(c.to_string());
            }

            out.extend(immediate);
        }

        out
    }

}
