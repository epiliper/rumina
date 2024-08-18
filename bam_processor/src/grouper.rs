extern crate bam;
use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::sync::Arc;

use std::str;
use strsim::hamming;

// this is the struct that contains functions used to group umis per the directional method
pub struct Grouper<'b> {
    pub umis: &'b Vec<String>,
}

impl<'b> Grouper<'b> {
    // gets the group neighbors of all group neighbors of a given UMI.
    pub fn depth_first_search(
        mut node: &'b String,
        adj_list: &IndexMap<&'b String, Vec<&'b String>>,
    ) -> HashSet<&'b String> {
        let mut searched: HashSet<&String> = HashSet::new();
        let mut queue: VecDeque<&String> = VecDeque::new();

        queue.push_back(node);
        searched.insert(node);

        while queue.len() > 0 {
            node = queue.pop_front().unwrap();
            if adj_list.contains_key(node) {
                &adj_list[node].iter().for_each(|next_node| {
                    if !searched.contains(next_node) {
                        queue.push_back(next_node);
                        searched.insert(next_node);
                    }
                });
            }
        }
        return searched;
    }

    // gets edit distance (hamming distance) between two umis
    pub fn edit_distance(ua: &String, ub: &String) -> usize {
        return hamming(ua, ub).unwrap();
    }

    pub fn get_substring_map(&self) -> IndexMap<&str, Vec<&String>> {
        let umi_length = self.umis.first().unwrap().len();
        let mid = umi_length / 2;

        let mut substring_map: IndexMap<&str, Vec<&String>> = IndexMap::new();

        self.umis.iter().for_each(|x| {
            let split = x.split_at(mid);
            substring_map.entry(split.0).or_insert(Vec::new()).push(x);
            substring_map.entry(split.1).or_insert(Vec::new()).push(x);
        });
        return substring_map;
    }

    pub fn iter_substring_neighbors(
        &self,
        substring_map: IndexMap<&'b str, Vec<&'b String>>,
    ) -> IndexMap<&'b String, IndexSet<&'b String>> {
        let mut neighbors: IndexMap<&'b String, IndexSet<&'b String>> = IndexMap::new();

        substring_map.iter().for_each(|x| {
            for umi in x.1 {
                neighbors.entry(umi).or_insert(IndexSet::new());
                let mut observed: IndexSet<&'b String> = IndexSet::new();

                let sub2 = (
                    substring_map.get(umi.split(x.0).next().unwrap()),
                    substring_map.get(umi.split(x.0).last().unwrap()),
                );
                match sub2 {
                    (Some(x), Some(y)) => {
                        observed.extend(x);
                        observed.extend(y);
                    }

                    (None, Some(y)) => {
                        observed.extend(y);
                    }

                    (Some(x), None) => {
                        observed.extend(x);
                    }

                    (None, None) => {}
                }

                neighbors[umi].extend(observed);
            }
        });

        return neighbors;
    }

    // groups umis via directional algorithm
    pub fn get_adj_list_directional(
        &self,
        counts: &HashMap<&String, i32>,
        substring_neighbors: IndexMap<&'b String, IndexSet<&'b String>>,
        threshold: usize,
    ) -> IndexMap<&'b String, Vec<&'b String>> {
        let mut adj_list: IndexMap<&'b String, Vec<&'b String>> =
            IndexMap::with_capacity(substring_neighbors.values().len());

        substring_neighbors.iter().for_each(|x| {
            let umi = x.0;
            let neighbors = x.1;
            adj_list.entry(umi).or_insert(Vec::new());

            for neighbor in neighbors {
                adj_list.entry(neighbor).or_insert(Vec::new());
                if Grouper::edit_distance(umi, neighbor) <= threshold && umi != neighbor {
                    if *counts.get(umi).unwrap() >= (counts.get(neighbor).unwrap() * 2 - 1) {
                        adj_list[umi].push(neighbor);
                    } else if *counts.get(neighbor).unwrap() >= (counts.get(umi).unwrap() * 2 - 1) {
                        adj_list[neighbor].push(umi);
                    }
                } else {
                }
            }
        });

        return adj_list;
    }

    // groups umis via bidirectional algorithm
    pub fn get_adj_list_bidirectional(
        &self,
        counts: &HashMap<&String, i32>,
        substring_neighbors: IndexMap<&'b String, IndexSet<&'b String>>,
        threshold: usize,
    ) -> IndexMap<&'b String, Vec<&'b String>> {
        let mut adj_list: IndexMap<&'b String, Vec<&'b String>> =
            IndexMap::with_capacity(substring_neighbors.values().len());

        // if a barcode is already part of a tree, don't group it again
        let mut found: HashSet<&String> = HashSet::with_capacity(adj_list.len());

        substring_neighbors.iter().for_each(|x| {
            let umi = x.0;
            let neighbors = x.1;

            if !found.contains(umi) {
                adj_list.entry(umi).or_insert(Vec::new());

                for neighbor in neighbors {
                    if !found.contains(neighbor) {
                        if Grouper::edit_distance(umi, neighbor) <= threshold && umi != neighbor {
                            if *counts.get(umi).unwrap() >= (counts.get(neighbor).unwrap() * 2 - 1)
                            {
                                adj_list[umi].push(neighbor);
                                found.insert(neighbor);
                            }
                        }
                    }
                }
            }
        });

        return adj_list;
    }

    // return a list of lists, comprising a UMI
    // with a list of grouped UMIs.
    // via depth-first-search
    // this is fed directly into the main_grouper function
    pub fn get_connected_components_par(
        &self,
        adj_list: IndexMap<&'b String, Vec<&'b String>>,
    ) -> Option<Vec<HashSet<&String>>> {
        let mut components: Vec<HashSet<&String>> = Vec::new();
        let mut found: HashSet<&String> = HashSet::new();

        if adj_list.len() > 0 {
            adj_list.keys().for_each(|node| {
                if !found.contains(node) {
                    let component = Grouper::depth_first_search(node, &adj_list);
                    found.extend(&component);
                    components.push(component);
                }
            });
            return Some(components);
        } else {
            return None;
        }
    }

    // get a list of UMIs, each with their own list of UMIs belonging to their group
    pub fn get_umi_groups(&self, clusters: Vec<HashSet<&'b String>>) -> Vec<Vec<&'b String>> {
        let mut observed: HashSet<&String> = HashSet::new();
        let mut groups: Vec<Vec<&String>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                observed.extend(&cluster);
                groups.push(cluster.into_iter().collect::<Vec<&String>>());
            } else {
                let mut temp_cluster: Vec<&String> = Vec::new();

                for node in cluster {
                    if !observed.contains(&node) {
                        temp_cluster.push(node);
                        observed.insert(node);
                    }
                }
                groups.push(temp_cluster);
            }
        }
        return groups;
    }

    pub fn no_clustering(
        &self,
        counts: HashMap<&'b String, i32>,
    ) -> (HashMap<&String, i32>, Option<Vec<Vec<&String>>>) {
        let umis = self
            .umis
            .iter()
            .map(|x| HashSet::from([x]))
            .collect::<Vec<HashSet<&'b String>>>();
        let final_umis = Some(self.get_umi_groups(umis));

        return (counts, final_umis);
    }

    pub fn cluster(
        &self,
        counts: HashMap<&'b String, i32>,
        grouping_method: Arc<&GroupingMethod>,
        num_umis: i64,
    ) -> (HashMap<&String, i32>, Option<Vec<Vec<&String>>>) {
        let clusterer = match *grouping_method {
            GroupingMethod::Directional => Grouper::get_adj_list_directional,
            GroupingMethod::Acyclic => Grouper::get_adj_list_bidirectional,
            GroupingMethod::Raw => {
                return self.no_clustering(counts);
            }
        };

        let substring_map = self.get_substring_map();
        let umis_to_compare = self.iter_substring_neighbors(substring_map);
        let adj_list = clusterer(&self, &counts, umis_to_compare, 1);
        let final_umis;

        if !adj_list.is_empty() {
            let clusters = self.get_connected_components_par(adj_list).unwrap();
            final_umis = Some(self.get_umi_groups(clusters));
        } else {
            final_umis = None;
        }
        return (counts, final_umis);
    }
}
