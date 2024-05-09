extern crate bam;
use indexmap::IndexMap;
use parking_lot::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::str;
use std::sync::Arc;
use strsim::hamming;

// this is the struct that contains functions used to group umis per the directional method
pub struct Processor<'b> {
    pub umis: &'b Vec<String>,
}

impl<'b> Processor<'b> {
    // for getting all unique umi keys from adjacency list
    // visits a UMI, gets all UMIs grouped with it, then moves onto next UMI.
    // doesn't add duplicate UMIs
    pub fn depth_first_search(
        mut node: &'b String,
        adj_list: &IndexMap<&'b String, HashSet<&'b String>>,
    ) -> VecDeque<&'b String> {
        let mut searched: VecDeque<&String> = VecDeque::new();
        let mut queue: VecDeque<&String> = VecDeque::new();

        queue.push_back(node);
        searched.push_back(node);

        while queue.len() > 0 {
            node = &queue.pop_front().unwrap();
            if adj_list.contains_key(node) {
                for next_node in &adj_list[node] {
                    if !searched.contains(&next_node) {
                        queue.push_back(next_node);
                        searched.push_back(next_node);
                    }
                }
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

        let substring_map: Arc<Mutex<IndexMap<&str, Vec<&String>>>> =
            Arc::new(Mutex::new(IndexMap::new()));

        self.umis.par_iter().for_each(|x| {
            let split = x.split_at(mid);
            substring_map
                .lock()
                .entry(split.0)
                .or_insert(Vec::new())
                .push(x);
            substring_map
                .lock()
                .entry(split.1)
                .or_insert(Vec::new())
                .push(x);
        });
        return Arc::try_unwrap(substring_map).unwrap().into_inner();
    }

    pub fn iter_substring_neighbors(
        &self,
        substring_map: IndexMap<&'b str, Vec<&'b String>>,
    ) -> IndexMap<&'b String, HashSet<&'b String>> {
        let neighbors: Arc<Mutex<IndexMap<&'b String, HashSet<&'b String>>>> =
            Arc::new(Mutex::new(IndexMap::new()));

        substring_map.par_iter().for_each(|x| {
            for umi in x.1 {
                neighbors.lock().entry(umi).or_insert(HashSet::new());
                let mut observed: HashSet<&'b String> = HashSet::new();

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

                neighbors.lock()[umi].extend(observed);
            }
        });

        return Arc::try_unwrap(neighbors).unwrap().into_inner();
    }

    // groups umis via directional algorithm

    pub fn get_adj_list_substring(
        &self,
        counts: HashMap<&String, i32>,
        substring_neighbors: IndexMap<&'b String, HashSet<&'b String>>,
        threshold: usize,
    ) -> IndexMap<&'b String, HashSet<&'b String>> {
        let adj_list: Arc<Mutex<IndexMap<&'b String, HashSet<&'b String>>>> = Arc::new(Mutex::new(
            IndexMap::with_capacity(substring_neighbors.values().len()),
        ));

        substring_neighbors.par_iter().for_each(|x| {
            let umi = x.0;
            let neighbors = x.1;
            adj_list.lock().entry(x.0).or_insert(HashSet::new());

            for neighbor in neighbors {
                adj_list.lock().entry(neighbor).or_insert(HashSet::new());
                if Processor::edit_distance(umi, neighbor) <= threshold && umi != neighbor {
                    if *counts.get(umi).unwrap() >= (counts.get(neighbor).unwrap() * 2 - 1) {
                        adj_list.lock()[umi].insert(neighbor);
                    } 
                    else if *counts.get(neighbor).unwrap() >= (counts.get(umi).unwrap() * 2 - 1) {
                        adj_list.lock()[umi].insert(neighbor);
                    }
                } else {
                }
            }
        });

        return Arc::try_unwrap(adj_list).unwrap().into_inner();
    }

    pub fn get_adj_list_substring_remove_singles(
        &self,
        counts: HashMap<&String, i32>,
        substring_neighbors: IndexMap<&'b String, HashSet<&'b String>>,
        threshold: usize,
    ) -> IndexMap<&'b String, HashSet<&'b String>> {
        let adj_list: Arc<Mutex<IndexMap<&'b String, HashSet<&'b String>>>> = Arc::new(Mutex::new(
            IndexMap::with_capacity(substring_neighbors.values().len()),
        ));

        substring_neighbors.par_iter().for_each(|x| {
            let umi = x.0;
            let neighbors = x.1;
            let umi_count = counts.get(umi).unwrap();
            if *umi_count == 1 {
            } else {

            adj_list.lock().entry(x.0).or_insert(HashSet::new());

            for neighbor in neighbors {
                // adj_list.lock().entry(neighbor).or_insert(HashSet::new());

                if Processor::edit_distance(umi, neighbor) <= threshold && umi != neighbor {
                    let neighbor_count = counts.get(neighbor).unwrap();
                    if  *umi_count >= (neighbor_count * 2 - 1) {
                        adj_list.lock()[umi].insert(neighbor);
                    } 
                    else if *neighbor_count > 1 && *neighbor_count >= (umi_count * 2 -1) {
                        adj_list.lock()[umi].insert(neighbor);
                    }
                } else {
                }
            }
            }
        });

        return Arc::try_unwrap(adj_list).unwrap().into_inner();
    }

    // return a list of lists, comprising a UMI
    // with a list of grouped UMIs.
    // via breadth-first-search
    // this is fed directly into the main_grouper function
    pub fn get_connected_components(
        &self,
        adj_list: IndexMap<&'b String, HashSet<&'b String>>,
    ) -> Option<Vec<VecDeque<&String>>> {
        let mut components: Vec<VecDeque<&String>> = Vec::new();
        let mut found: Vec<&String> = Vec::new();

        if adj_list.len() > 0 {
            for node in adj_list.keys() {
                if !found.contains(node) {
                    let component = Processor::depth_first_search(node, &adj_list);
                    found.extend(&component);
                    components.push(component);
                }
            }
            return Some(components);
        } else {
            return None;
        }
    }

    pub fn get_connected_components_par(
        &self,
        adj_list: IndexMap<&'b String, HashSet<&'b String>>,
    ) -> Option<Vec<VecDeque<&String>>> {
        let mut components: Arc<Mutex<Vec<VecDeque<&String>>>> = Arc::new(Mutex::new(Vec::new()));
        let mut found: Arc<Mutex<HashSet<&String>>> = Arc::new(Mutex::new(HashSet::new()));

        if adj_list.len() > 0 {
            adj_list.par_keys().for_each(
                |node| {
                if !found.lock().contains(node) {
                    let component = Processor::depth_first_search(node, &adj_list);
                    found.lock().extend(&component);
                    components.lock().push(component);
                }
                }
            );
            return Some(Arc::try_unwrap(components).unwrap().into_inner());
        } else {
            return None;
        }
    }



    // get a list of UMIs, each with their own list of UMIs belonging to their group
    pub fn group_directional(&self, clusters: Vec<VecDeque<&'b String>>) -> Vec<Vec<&'b String>> {
        let mut observed: Vec<&String> = Vec::new();
        let mut groups: Vec<Vec<&String>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                observed.push(&cluster.get(0).unwrap());
                groups.push(cluster.into());
            } else {
                let mut temp_cluster: Vec<&String> = Vec::new();

                for node in cluster {
                    if !observed.contains(&&node) {
                        temp_cluster.push(node);
                        observed.push(node);
                    }
                }
                groups.push(temp_cluster);
            }
        }
        return groups;
    }

    // driver code for directional method,
    // and UMI organization and grouping
    pub fn main_grouper(&self, counts: HashMap<&String, i32>) -> Option<Vec<Vec<&String>>> {
        let substring_map = self.get_substring_map();
        let neighbors = self.iter_substring_neighbors(substring_map);
        // let directional_output = self.get_adj_list_substring(counts, neighbors, 1);
        let directional_output = self.get_adj_list_substring_remove_singles(counts, neighbors, 1);
        let adj_list = directional_output;
        let final_umis;
        println!{"len adj_list {}", adj_list.len()}
        if adj_list.len() > 0 {
            let clusters = self.get_connected_components_par(adj_list).unwrap();
            final_umis = Some(self.group_directional(clusters));
        } else {
            final_umis = None;
        }

        return final_umis;
    }
}
