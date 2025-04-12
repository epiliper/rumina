use crate::group::bktree::NGramBKTree;
use crate::ngram::ngram;
use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use std::cell::{Ref, RefCell};
use std::collections::HashMap;
use std::sync::Arc;

pub struct Grouper<'a> {
    pub umis: &'a Vec<String>,
    pub max_edit: usize,
    pub umi_len: usize,
    pub ngram_maker: ngram::NgramMaker<'a>,
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

    pub fn ngrams(&self, umi: &'a str) -> Ref<Vec<&str>> {
        self.ngram_maker.ngrams(umi)
    }

    pub fn get_substring_map(&self) -> IndexMap<&str, IndexSet<&str>> {
        let mut substring_map: IndexMap<&str, IndexSet<&str>> = IndexMap::new();

        for umi in self.umis.iter() {
            self.ngrams(umi).iter().for_each(|slice| {
                substring_map
                    .entry(&slice)
                    .or_insert(IndexSet::new())
                    .insert(umi);
            })
        }

        substring_map
    }

    pub fn cluster(
        &'a self,
        counts: HashMap<&'a str, i32>,
        grouping_method: Arc<&GroupingMethod>,
    ) -> Option<Vec<IndexSet<String>>> {
        let substring_neighbors = self.get_substring_map();

        let mut out = Vec::new();

        let mut bktree = NGramBKTree::init_empty(None, self.umi_len);
        bktree.count_map = RefCell::new(counts.clone());
        for (ng, nei) in substring_neighbors {
            bktree.populate_from_neighbors(ng, nei);
        }

        self.umis.iter().for_each(|u| {
            let max_count = counts.get(u.as_str()).unwrap() / 2 + 1;
            let o = bktree.remove_near(u.as_str(), 1, max_count, &self.ngram_maker);
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
