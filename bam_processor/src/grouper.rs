use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::iter::zip;
use std::sync::Arc;

use std::str;
use strsim::hamming;

// this is the struct that contains functions used to group umis per the directional method
pub struct Grouper<'b> {
    pub umis: &'b Vec<String>,
}

impl<'b> Grouper<'b> {
    // gets the group neighbors of all group neighbors of a given UMI.
    pub fn breadth_first_search(
        mut node: &'b String,
        adj_list: &IndexMap<&'b String, Vec<&'b String>>,
    ) -> HashSet<&'b String> {
        let mut searched: HashSet<&String> = HashSet::new();
        let mut queue: VecDeque<&String> = VecDeque::new();

        queue.push_back(node);
        searched.insert(node);

        while !queue.is_empty() {
            node = queue.pop_front().unwrap();
            if adj_list.contains_key(node) {
                adj_list[node].iter().for_each(|next_node| {
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

    pub fn get_substring_map(&self) -> IndexMap<(usize, usize), IndexMap<&str, Vec<&String>>> {
        let umi_length = self.umis.first().unwrap().len();
        let mut slices = Vec::with_capacity(umi_length * 2);

        let mut offset = 0;
        let chunk_size = umi_length / 2;
        let mut remainder = umi_length % 2;

        for _ in 0..2 {
            let end = offset + chunk_size + (remainder > 0) as usize;
            slices.push((offset, end));
            offset = end;
            if remainder > 0 {
                remainder -= 1;
            }
        }

        let mut substring_map: IndexMap<(usize, usize), IndexMap<&str, Vec<&String>>> =
            IndexMap::new();

        let indices = slices.clone();

        for (slice, idx) in zip(slices, indices) {
            let slice_entry = substring_map.entry(idx).or_insert_with(|| IndexMap::new());

            for umi in self.umis {
                slice_entry
                    .entry(&umi[slice.0..slice.1])
                    .or_insert_with(|| vec![umi])
                    .push(umi);
            }
        }

        return substring_map;
    }

    pub fn iter_substring_neighbors(
        &self,
        substring_map: IndexMap<(usize, usize), IndexMap<&'b str, Vec<&'b String>>>,
    ) -> IndexMap<&'b String, IndexSet<&'b String>> {
        let mut neighbors: IndexMap<&'b String, IndexSet<&'b String>> = IndexMap::new();
        for (i, u) in self.umis.iter().enumerate() {
            let mut observed: HashSet<&String> = HashSet::new();

            neighbors.entry(u).or_insert(IndexSet::new());

            for (slice, substrings) in &substring_map {
                neighbors[u].extend(substrings.get(&u[slice.0..slice.1]).unwrap())
            }

            neighbors[u].retain(|nbr| !observed.contains(nbr));
            observed.insert(u);
        }

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
    pub fn get_connected_components(
        &self,
        adj_list: IndexMap<&'b String, Vec<&'b String>>,
    ) -> Option<Vec<HashSet<&String>>> {
        let mut components: Vec<HashSet<&String>> = Vec::new();
        let mut found: HashSet<&String> = HashSet::new();

        if adj_list.len() > 0 {
            adj_list.keys().for_each(|node| {
                if !found.contains(node) {
                    let component = Grouper::breadth_first_search(node, &adj_list);
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
            let mut temp_cluster: Vec<&String> = Vec::new();

            for node in cluster {
                if !observed.contains(&node) {
                    temp_cluster.push(node);
                    observed.insert(node);
                }
            }
            groups.push(temp_cluster);
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
            let clusters = self.get_connected_components(adj_list).unwrap();
            final_umis = Some(self.get_umi_groups(clusters));
        } else {
            final_umis = None;
        }
        return (counts, final_umis);
    }
}
