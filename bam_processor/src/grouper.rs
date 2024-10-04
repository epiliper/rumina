use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::iter::zip;
use std::sync::Arc;

use std::str;
use strsim::hamming;

// gets edit distance (hamming distance) between two umis
pub fn edit_distance(ua: &str, ub: &str) -> usize {
    hamming(ua, ub).unwrap()
}

// this struct contains logic for error correction of UMIs downstream of read batching.
// this phase of the pipeline is driven by the cluster() and no_cluster() functions.
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
                adj_list.get(node).unwrap().iter().for_each(|next_node| {
                    if !searched.contains(next_node) {
                        queue.push_back(next_node);
                        searched.insert(next_node);
                    }
                });
            }
        }
        searched
    }

    // groups UMIs with substring neighbors to reduce number of comparisons.
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
            let slice_entry = substring_map.entry(idx).or_default();

            for umi in self.umis {
                slice_entry
                    .entry(&umi[slice.0..slice.1])
                    .or_insert(vec![umi])
                    .push(umi);
            }
        }

        substring_map
    }

    // create an IndexMap (k, [v]) where UMI k has edit distance calculated for each UMI within [v].
    pub fn iter_substring_neighbors(
        &self,
        substring_map: IndexMap<(usize, usize), IndexMap<&'b str, Vec<&'b String>>>,
    ) -> IndexMap<&'b String, IndexSet<&'b String>> {
        let mut neighbors: IndexMap<&'b String, IndexSet<&'b String>> = IndexMap::new();

        let mut observed: HashSet<&String> = HashSet::new();
        for u in self.umis.iter() {
            for (slice, substrings) in &substring_map {
                neighbors
                    .entry(u)
                    .or_insert_with(|| IndexSet::new())
                    .extend(substrings.get(&u[slice.0..slice.1]).unwrap())
            }

            observed.insert(u);
            neighbors
                .get_mut(u)
                .unwrap()
                .retain(|nbr| !observed.contains(nbr));
        }

        neighbors
    }

    // if UMI group A within edit distance threshold and has ≥2n-1 read count compared to group B,
    // then (group A) ---> (group B)
    pub fn add_edge_directional(
        &self,
        umi_a: &String,
        umi_b: &String,
        counts: &HashMap<&String, i32>,
        threshold: usize,
    ) -> bool {
        edit_distance(umi_a, umi_b) <= threshold
            && umi_a != umi_b
            && *counts.get(umi_a).unwrap() >= counts.get(umi_b).unwrap() * 2 - 1
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

            adj_list.entry(umi).or_default();

            for neighbor in neighbors {
                if self.add_edge_directional(umi, neighbor, counts, threshold) {
                    adj_list.entry(umi).or_insert(vec![neighbor]).push(neighbor);
                }

                if self.add_edge_directional(neighbor, umi, counts, threshold) {
                    adj_list.entry(neighbor).or_insert(vec![umi]).push(umi);
                }
            }
        });
        adj_list
    }

    // groups umis via acyclic algorithm
    pub fn get_adj_list_acyclic(
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
                adj_list.entry(umi).or_default();

                for neighbor in neighbors {
                    if !found.contains(neighbor)
                        && self.add_edge_directional(umi, neighbor, counts, threshold)
                    {
                        adj_list[umi].push(neighbor);
                        found.insert(neighbor);
                    }
                }
            }
        });
        adj_list
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

        if !adj_list.is_empty() {
            adj_list.keys().for_each(|node| {
                if !found.contains(node) {
                    let component = Grouper::breadth_first_search(node, &adj_list);
                    found.extend(&component);
                    components.push(component);
                }
            });

            Some(components)
        } else {
            None
        }
    }

    // get a list of UMIs, each with their own list of UMIs belonging to their group
    pub fn get_umi_groups(&self, clusters: Vec<HashSet<&'b String>>) -> Vec<Vec<&'b String>> {
        let mut observed: HashSet<&String> = HashSet::new();
        let mut groups: Vec<Vec<&String>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                let node = cluster.iter().next().unwrap();
                groups.push(vec![node]);
                observed.insert(node);
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
        groups
    }

    // used with the raw method. No UMI error correction.
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

        (counts, final_umis)
    }

    // used with directional and acyclic methods. Driver function of UMI-error correction
    pub fn cluster(
        &self,
        counts: HashMap<&'b String, i32>,
        grouping_method: Arc<&GroupingMethod>,
    ) -> (HashMap<&String, i32>, Option<Vec<Vec<&String>>>) {
        let clusterer = match *grouping_method {
            GroupingMethod::Directional => Grouper::get_adj_list_directional,
            GroupingMethod::Acyclic => Grouper::get_adj_list_acyclic,
            GroupingMethod::Raw => {
                return self.no_clustering(counts);
            }
        };

        let substring_map = self.get_substring_map();
        let umis_to_compare = self.iter_substring_neighbors(substring_map);
        let adj_list = clusterer(self, &counts, umis_to_compare, 1);

        if !adj_list.is_empty() {
            let clusters = self.get_connected_components(adj_list).unwrap();
            let final_umis = Some(self.get_umi_groups(clusters));
            (counts, final_umis)
        } else {
            (counts, None)
        }
    }
}
