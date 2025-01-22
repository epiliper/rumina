use crate::GroupingMethod;
use indexmap::{IndexMap, IndexSet};
use log::{debug, warn};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::iter::zip;
use std::sync::Arc;
use std::vec::IntoIter;

use std::str;
use strsim::hamming;

// gets edit distance (hamming distance) between two umis
pub fn edit_distance(ua: &str, ub: &str) -> usize {
    hamming(ua, ub).unwrap()
}

#[derive(Debug)]
pub struct GroupIterator<'b> {
    pub clusters: IntoIter<HashSet<&'b str>>,
    observed: HashSet<&'b str>,
}

impl<'b> Iterator for GroupIterator<'b> {
    type Item = Vec<&'b str>;

    fn next(&mut self) -> Option<Self::Item> {
        let next = self.clusters.next();

        if next.is_some() {
            let cluster = next.unwrap();
            if cluster.len() == 1 {
                let node = cluster.iter().next().unwrap();
                self.observed.insert(node);
                return Some(vec![node]);
            }

            let mut temp_cluster: Vec<&str> = Vec::new();

            for node in cluster {
                if !self.observed.contains(&node) {
                    temp_cluster.push(node);
                    self.observed.insert(node);
                }
            }

            return Some(temp_cluster);
        }

        None
    }
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
            warn!("Empty adjacency list");
            None
        }
    }

    // used with the raw method. No UMI error correction.
    pub fn no_clustering(
        &self,
        counts: HashMap<&'b str, i32>,
    ) -> (HashMap<&str, i32>, GroupIterator) {
        let umis = self
            .umis
            .iter()
            .map(|x| HashSet::from([x.as_str()]))
            .collect::<Vec<HashSet<&'b str>>>();
        // let final_umis = self.get_umi_groups(umis);
        let final_umis = GroupIterator {
            clusters: umis.into_iter(),
            observed: HashSet::new(),
        };

        (counts, final_umis)
    }

    // used with directional and acyclic methods. Driver function of UMI-error correction
    pub fn cluster(
        &self,
        counts: HashMap<&'b str, i32>,
        grouping_method: Arc<&GroupingMethod>,
    ) -> (HashMap<&str, i32>, Option<GroupIterator>) {
        let clusterer = match *grouping_method {
            GroupingMethod::Directional => Grouper::get_adj_list_directional,
            GroupingMethod::Acyclic => Grouper::get_adj_list_acyclic,
            GroupingMethod::Raw => {
                let (counts, groups) = self.no_clustering(counts);
                return (counts, Some(groups));
            }
        };

        let substring_map = self.get_substring_map();
        let umis_to_compare = self.iter_substring_neighbors(substring_map);
        let adj_list = clusterer(self, &counts, umis_to_compare, 1);
        debug!(
            "Generated adjacency list with {} entries...",
            adj_list.len()
        );

        if !adj_list.is_empty() {
            let clusters = self.get_connected_components(adj_list).unwrap();
            // let final_umis = self.get_umi_groups(clusters);
            let final_umis = GroupIterator {
                clusters: clusters.into_iter(),
                observed: HashSet::new(),
            };
            (counts, Some(final_umis))
        } else {
            (counts, None)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_edit_distance() {
        let ua = "ATCG";
        let ub = "ATGG";
        assert_eq!(edit_distance(ua, ub), 1);

        let ua = "ATCG";
        let ub = "ATCG";
        assert_eq!(edit_distance(ua, ub), 0);

        let ua = "ATCG";
        let ub = "GGGG";
        assert_eq!(edit_distance(ua, ub), 3);
    }

    #[test]
    fn test_get_substring_map() {
        let umis = vec!["ATCG".to_string(), "ATGG".to_string(), "ACGG".to_string()];
        let grouper = Grouper { umis: &umis };
        let substring_map = grouper.get_substring_map();

        // Check if the substring map is correctly generated
        assert!(!substring_map.is_empty());
        for (_slice, substrings) in &substring_map {
            for (substring, umis) in substrings {
                assert!(!umis.is_empty());
                for umi in umis {
                    assert!(umi.contains(substring));
                }
            }
        }

        assert_eq!(
            *substring_map.get(&(0, 2)).unwrap().get("AT").unwrap(),
            vec!["ATCG", "ATGG"]
        )
    }

    #[test]
    fn test_add_edge_directional() {
        let umis = vec!["ATCG".to_string(), "ATGG".to_string(), "ACGG".to_string()];
        let grouper = Grouper { umis: &umis };
        let mut counts = HashMap::new();
        counts.insert("ATCG", 10);
        counts.insert("ATGG", 5);
        counts.insert("ACGG", 3);

        assert!(grouper.add_edge_directional("ATCG", "ATGG", &counts, 1));
        assert!(!grouper.add_edge_directional("ATGG", "ATCG", &counts, 1));
        assert!(!grouper.add_edge_directional("ATCG", "ACGG", &counts, 1));
    }

    #[test]
    fn test_get_adj_list_directional() {
        let umis = vec!["ATCG".to_string(), "ATGG".to_string(), "ACGG".to_string()];
        let grouper = Grouper { umis: &umis };
        let mut counts = HashMap::new();
        counts.insert("ATCG", 10);
        counts.insert("ATGG", 5);
        counts.insert("ACGG", 3);

        let substring_map = grouper.get_substring_map();
        let neighbors = grouper.iter_substring_neighbors(substring_map);
        let adj_list = grouper.get_adj_list_directional(&counts, neighbors, 1);

        // Check if the adjacency list is correctly generated
        assert!(!adj_list.is_empty());
        assert_eq!(adj_list.get("ATCG").unwrap()[0], "ATGG");
        assert_eq!(adj_list.get("ATGG").unwrap()[0], "ACGG");
    }

    #[test]
    fn test_get_umi_groups() {
        let umis = vec![
            "ATCG".to_string(),
            "ATGG".to_string(),
            "ACGG".to_string(),
            "GATT".to_string(),
        ];
        let grouper = Grouper { umis: &umis };
        let mut counts = HashMap::new();
        counts.insert("ATCG", 10);
        counts.insert("ATGG", 5);
        counts.insert("ACGG", 3);
        counts.insert("GATT", 25);

        let substring_map = grouper.get_substring_map();
        let neighbors = grouper.iter_substring_neighbors(substring_map);
        let adj_list = grouper.get_adj_list_directional(&counts, neighbors, 1);
        let components = grouper.get_connected_components(adj_list).unwrap();
        // let groups = grouper.get_umi_groups(components);
        let groups = GroupIterator {
            clusters: components.into_iter(),
            observed: HashSet::new(),
        };

        let groups = groups.into_iter().collect::<Vec<_>>();

        // Check if UMI groups are correctly generated
        println!("{:?}", groups);
        assert!(!groups.is_empty());
        assert_eq!(groups.len(), 2);
        for group in &groups {
            assert!(!group.is_empty());
        }
    }

    #[test]
    fn test_cluster_directional() {
        let umis = vec![
            "ATCG".to_string(),
            "ATGG".to_string(),
            "ACGG".to_string(),
            "GATT".to_string(),
        ];
        let grouper = Grouper { umis: &umis };
        let mut counts = HashMap::new();
        counts.insert("ATCG", 10);
        counts.insert("ATGG", 5);
        counts.insert("ACGG", 3);
        counts.insert("GATT", 25);

        let grouping_method = Arc::new(&GroupingMethod::Directional);
        let (_counts, groups) = grouper.cluster(counts, grouping_method);
        let groups = groups.unwrap().into_iter().collect::<Vec<_>>();
        println!("{:?}", groups);

        // Check if clustering is correctly performed
        let groups = groups;
        assert_eq!(groups.len(), 2)
    }

    #[test]
    fn test_no_clustering() {
        let umis = vec!["ATCG".to_string(), "ATGG".to_string(), "ACGG".to_string()];
        let grouper = Grouper { umis: &umis };
        let mut counts = HashMap::new();
        counts.insert("ATCG", 10);
        counts.insert("ATGG", 5);
        counts.insert("ACGG", 3);

        let (_counts, groups) = grouper.no_clustering(counts);
        let groups = groups.collect::<Vec<_>>();

        // Check if no clustering is correctly performed
        assert_eq!(groups.len(), umis.len());
        for group in groups {
            assert_eq!(group.len(), 1);
        }
    }
}
