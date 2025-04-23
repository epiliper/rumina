use crate::group::bktree::NGramBKTree;
use crate::ngram::NgramMaker;
use crate::processor::UmiHistogram;
use crate::GroupingMethod;
use indexmap::IndexSet;
use smol_str::SmolStr;
use std::collections::VecDeque;
use std::sync::Arc;

pub struct Grouper<'a> {
    pub umis: &'a Vec<SmolStr>,
    pub ngram_maker: NgramMaker,
    percentage: f32,
    max_edit: u32,
}

impl<'a> Grouper<'a> {
    pub fn new(umis: &'a Vec<SmolStr>, max_edit: u32, percentage: f32, umi_len: usize) -> Self {
        assert!(percentage > 0.0 && percentage <= 1.0);
        let ngram_maker = NgramMaker::new(
            (max_edit + 1)
                .try_into()
                .expect("Unable to convert value for ngram_maker"),
            umi_len,
        );

        Self {
            umis,
            ngram_maker,
            percentage,
            max_edit,
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
        mut counts: UmiHistogram<'a>,
        grouping_method: Arc<&'a GroupingMethod>,
    ) -> impl Iterator<Item = IndexSet<SmolStr>> + 'a {
        let mut bk = match *grouping_method {
            GroupingMethod::Directional | GroupingMethod::Acyclic => {
                Some(self.init_bktree(&counts))
            }
            GroupingMethod::Raw => None,
        };

        let process = move |u: &'a SmolStr| -> IndexSet<SmolStr> {
            match *grouping_method {
                GroupingMethod::Directional => {
                    self.visit_and_remove_all(u, self.max_edit, &mut counts, bk.as_mut().unwrap())
                }
                GroupingMethod::Acyclic => {
                    self.visit_and_remove_immediate(u, 1, &mut counts, bk.as_mut().unwrap())
                }
                GroupingMethod::Raw => self.remove_single(u),
            }
        };

        self.umis.iter().map(process).filter(|o| !o.is_empty())
    }

    pub fn init_bktree(&self, counts: &UmiHistogram) -> NGramBKTree {
        let mut bktree = NGramBKTree::init_empty(None);

        self.umis.iter().for_each(|u| {
            let count = *counts
                .get(u.as_str())
                .expect("FATAL ERROR: umi to deduplicate is not found in count map");
            bktree.populate_single(u, count.0, &self.ngram_maker);
        });
        bktree
    }

    pub fn remove_single(&self, umi: &str) -> IndexSet<SmolStr> {
        IndexSet::from([SmolStr::from(umi)])
    }

    pub fn visit_and_remove_immediate(
        &self,
        umi: &str,
        k: u32,
        counts: &mut UmiHistogram,
        bktree: &mut NGramBKTree,
    ) -> IndexSet<SmolStr> {
        let max_count =
            (self.percentage * (counts.get(umi).expect("Count map invalid").0 + 1) as f32) as i32;
        bktree.remove_near(umi, k, max_count, &self.ngram_maker, counts)
    }

    /// for every queried node, cluster its immediate offshoots distant by <= k edits and the given count percentage, then retrieve the
    /// offshoots of each offshoot, until no more offshoots are found.
    pub fn visit_and_remove_all(
        &self,
        umi: &str,
        k: u32,
        counts: &mut UmiHistogram,
        bktree: &mut NGramBKTree,
    ) -> IndexSet<SmolStr> {
        let mut to_cluster = VecDeque::from([umi.to_string()]);
        let mut out = IndexSet::new();

        while let Some(root) = to_cluster.pop_front() {
            // TODO: make percentage adjustable
            let max_count =
                (self.percentage * (counts.get(root.as_str()).unwrap().0 + 1) as f32) as i32;
            let immediate =
                bktree.remove_near(root.as_str(), k, max_count, &self.ngram_maker, counts);

            // if directional, we limit a umi's cluster to a depth of one maximum
            for c in immediate.iter().filter(|c| **c != root) {
                to_cluster.push_back(c.to_string());
            }

            out.extend(immediate);
        }

        out
    }
}
