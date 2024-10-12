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
pub struct Grouper<'a> {
    pub umis: &'a Vec<String>,
}

impl<'b> Grouper<'b> {
    // gets the group neighbors of all group neighbors of a given UMI.
    pub fn breadth_first_search(
        mut node: &'b str,
        adj_list: &IndexMap<&str, Vec<&'b str>>,
    ) -> HashSet<&'b str> {
        let mut searched: HashSet<&str> = HashSet::new();
        let mut queue: VecDeque<&str> = VecDeque::new();

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
    pub fn get_substring_map(&self) -> IndexMap<(usize, usize), IndexMap<&str, Vec<&str>>> {
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

        let mut substring_map: IndexMap<(usize, usize), IndexMap<&str, Vec<&str>>> =
            IndexMap::new();

        let indices = slices.clone();

        for (slice, idx) in zip(slices, indices) {
            let slice_entry = substring_map.entry(idx).or_default();

            for umi in self.umis.iter() {
                slice_entry
                    .entry(&umi[slice.0..slice.1])
                    .or_insert(Vec::new())
                    .push(umi);
            }
        }

        substring_map
    }

    // create an IndexMap (k, [v]) where UMI k has edit distance calculated for each UMI within [v].
    pub fn iter_substring_neighbors<'c>(
        &'c self,
        substring_map: IndexMap<(usize, usize), IndexMap<&'c str, Vec<&'c str>>>,
    ) -> impl Iterator<Item = (&str, IndexSet<&str>)> {
        self.umis.iter().map(move |u| {
            let mut neighbors: IndexSet<&str> = IndexSet::new();
            for (slice, substrings) in &substring_map {
                neighbors.extend(substrings.get(&u[slice.0..slice.1]).unwrap());
            }
            (u.as_str(), neighbors)
        })
    }

    // if UMI group A within edit distance threshold and has ≥2n-1 read count compared to group B,
    // then (group A) ---> (group B)
    pub fn add_edge_directional(
        &self,
        umi_a: &str,
        umi_b: &str,
        counts: &HashMap<&str, i32>,
        threshold: usize,
    ) -> bool {
        edit_distance(umi_a, umi_b) <= threshold
            && umi_a != umi_b
            && *counts.get(umi_a).unwrap() >= counts.get(umi_b).unwrap() * 2 - 1
    }

    // groups umis via directional algorithm
    pub fn get_adj_list_directional<'d>(
        &'d self,
        counts: &HashMap<&str, i32>,
        substring_neighbors: impl Iterator<Item = (&'d str, IndexSet<&'d str>)>,
        threshold: usize,
    ) -> IndexMap<&'d str, Vec<&'d str>> {
        let mut adj_list: IndexMap<&'d str, Vec<&'d str>> = IndexMap::new();

        substring_neighbors.for_each(|x| {
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
    pub fn get_adj_list_acyclic<'d>(
        &'d self,
        counts: &HashMap<&str, i32>,
        substring_neighbors: impl Iterator<Item = (&'d str, IndexSet<&'d str>)>,
        threshold: usize,
    ) -> IndexMap<&'d str, Vec<&'d str>> {
        let mut adj_list: IndexMap<&'d str, Vec<&'d str>> = IndexMap::new();

        // if a barcode is already part of a tree, don't group it again
        let mut found: HashSet<&str> = HashSet::new();

        substring_neighbors.for_each(|x| {
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
    pub fn get_connected_components<'e>(
        &'e self,
        adj_list: IndexMap<&'e str, Vec<&'e str>>,
    ) -> Option<Vec<HashSet<&str>>> {
        let mut components: Vec<HashSet<&str>> = Vec::new();
        let mut found: HashSet<&str> = HashSet::new();

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
    pub fn get_umi_groups(&self, clusters: Vec<HashSet<&'b str>>) -> Vec<Vec<&'b str>> {
        let mut observed: HashSet<&str> = HashSet::new();
        let mut groups: Vec<Vec<&str>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                let node = cluster.iter().next().unwrap();
                groups.push(vec![node]);
                observed.insert(node);
            } else {
                let mut temp_cluster: Vec<&str> = Vec::new();

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
        counts: HashMap<&'b str, i32>,
    ) -> (HashMap<&str, i32>, Option<Vec<Vec<&str>>>) {
        let umis = self
            .umis
            .iter()
            .map(|x| HashSet::from([x.as_str()]))
            .collect::<Vec<HashSet<&'b str>>>();
        let final_umis = Some(self.get_umi_groups(umis));

        (counts, final_umis)
    }

    // used with directional and acyclic methods. Driver function of UMI-error correction
    pub fn cluster(
        &self,
        counts: HashMap<&'b str, i32>,
        grouping_method: Arc<&GroupingMethod>,
    ) -> (HashMap<&str, i32>, Option<Vec<Vec<&str>>>) {
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
