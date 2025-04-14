use crate::group::bktree::NGramBKTree;
use crate::ngram::ngram;
use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use std::cell::{Ref, RefCell, RefMut};
use std::collections::{HashMap, VecDeque};
use std::sync::Arc;

pub struct Grouper<'a> {
    pub umis: &'a Vec<String>,
    pub max_edit: usize,
    pub umi_len: usize,
    pub ngram_maker: ngram::NgramMaker,
}

impl<'a> Grouper<'a> {
    pub fn new(umis: &'a Vec<String>, max_edit: usize, umi_len: usize) -> Self {
        let ngram_maker = ngram::NgramMaker::new(max_edit + 1, umi_len);

        Self {
            umis,
            max_edit,
            umi_len,
            ngram_maker,
        }
    }

    pub fn ngrams(&self, umi: &'a str) -> RefMut<Vec<String>> {
        self.ngram_maker.ngrams(umi)
    }

    pub fn get_substring_map(&self) -> IndexMap<String, IndexSet<&str>> {
        let mut substring_map: IndexMap<String, IndexSet<&str>> = IndexMap::new();

        for umi in self.umis.iter() {
            self.ngrams(umi).iter().for_each(|slice| {
                substring_map
                    .entry(slice.to_string())
                    .or_insert(IndexSet::new())
                    .insert(umi);
            })
        }

        substring_map
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
            let max_count = counts.get(root.as_str()).unwrap() / 2 + 1;
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
    ) -> Option<Vec<IndexSet<String>>> {
        let substring_neighbors = self.get_substring_map();

        let mut out = Vec::new();

        let mut bktree = NGramBKTree::init_empty(None, self.umi_len);
        bktree.count_map = counts.clone();
        let original_len = substring_neighbors.len();

        for (ng, nei) in substring_neighbors {
            bktree.populate_from_neighbors(ng, nei);
        }

        assert_eq!(bktree.ngram_tree_map.len(), original_len);

        self.umis.iter().for_each(|u| {
            let o = self.visit_and_remove_all(u, 1, &counts, &mut bktree);
            if !o.is_empty() {
                out.push(o);
            }
        });

        if out.is_empty() {
            None
        } else {
            Some(out)
        }
    }
}
