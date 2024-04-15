extern crate bam;
extern crate rust_htslib;
use indexmap::IndexMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

pub fn breadth_first_search(
    mut node: String,
    adj_list: &IndexMap<String, HashSet<String>>,
    // ) -> IndexMap<String, HashSet<String>>
) -> VecDeque<String> {
    // let mut searched: IndexMap<String, HashSet<String>> = IndexMap::new();
    // let mut queue: IndexMap<String, HashSet<String>> = IndexMap::new();
    let mut searched: VecDeque<String> = VecDeque::new();
    let mut queue: VecDeque<String> = VecDeque::new();

    queue.push_back(node);

    while queue.len() > 0 {
        node = queue.pop_front().unwrap();
        for next_node in &adj_list[&node] {
            if !searched.contains(next_node) {
                queue.push_back(next_node.to_string());
                searched.push_back(next_node.to_string());
            }
        }
    }
    return searched;
}

pub fn edit_distance(ua: &String, ub: &String) -> i32 {
    let la = ua.len();
    let lb = ub.len();
    let mut edit_distance: i32 = 0;

    if la != lb {
        return i32::MAX;
    } else {
        for i in 0..la {
            if ua.chars().collect::<Vec<_>>()[i] != ub.chars().collect::<Vec<_>>()[i] {
                edit_distance += 1;
            } else {
                continue;
            }
        }
        return edit_distance;
    }
}

pub struct NeighboringUMIs {
    pub neighbors: HashMap<String, String>,
}

impl NeighboringUMIs {
    pub fn iter(&self) -> NeighborIterator {
        NeighborIterator {
            nbrs: &self.neighbors,
            index: 0,
        }
    }
}

pub struct NeighborIterator<'a> {
    nbrs: &'a HashMap<String, String>,
    index: i32,
}

impl<'a> Iterator for NeighborIterator<'a> {
    type Item = (&'a String, &'a String);

    fn next(&mut self) -> Option<Self::Item> {
        while self.index < self.nbrs.len().try_into().unwrap() {
            for (substring, umi) in self.nbrs {
                self.index += 1;
                return Some((substring, umi));
            }
        }
        return None;
    }
}

pub struct Processor {
    // pub groups: UMIReads,
}
impl Processor {
    // create a list of indexes that split UMIs into substrings
    pub fn get_substr_slices(umi_length: i32, threshold: i32) -> Option<Vec<i32>> {
        let dividend = umi_length / threshold;
        let remainder = umi_length % threshold;
        let sub_sizes = vec![
            (dividend + 1) * remainder,
            (dividend * (threshold - remainder)),
        ];
        let mut offset = 0;
        let mut slices = Vec::new();
        if sub_sizes.is_empty() {
            return None;
        } else {
            for s in sub_sizes {
                slices.push(offset + s);
                offset += 1;
            }
            return Some(slices);
        }
    }

    // use indices to split UMIs, return a dict comparing all UMIs with similar substrings
    pub fn build_substr_idx(
        mut umis: Vec<String>,
        umi_length: i32,
        threshold: i32,
    ) -> HashMap<i32, HashMap<String, HashSet<String>>> {
        let mut substr_idx: HashMap<i32, HashMap<String, HashSet<String>>> = HashMap::new();
        let slices = Processor::get_substr_slices(umi_length, threshold + 1).unwrap();

        for idx in slices {
            for u in &mut umis {
                let u_sub = u.split_off(idx.try_into().unwrap());
                substr_idx
                    .entry(idx)
                    .or_default()
                    .entry(u_sub.to_string())
                    .or_default()
                    .insert(u.to_string());
            }
        }

        return substr_idx;
    }

    pub fn iter_nearest_neighbors(
        mut umis: Vec<String>,
        substr_idx: HashMap<i32, HashMap<String, HashSet<String>>>,
    ) -> Option<NeighboringUMIs> {
        let old = umis.clone();
        let mut neighboring_umis = NeighboringUMIs {
            neighbors: HashMap::new(),
        };
        if !umis.is_empty() {
            for (i, umi) in umis.iter_mut().enumerate() {
                for (idx, substring) in &substr_idx {
                    let usub = &umi.split_off((*idx).try_into().unwrap());
                    for u in &substring[usub] {
                        if !&neighboring_umis.neighbors.contains_key(u) & old[..i].contains(u) {
                            neighboring_umis
                                .neighbors
                                .insert(u.to_string(), umi.to_string());
                        } else {
                            continue;
                        }
                    }
                }
            }
            return Some(neighboring_umis);
        } else {
            return None;
        }
    }

    pub fn get_adj_list_directional(
        &self,
        umis: Vec<String>,
        counts: HashMap<String, i32>,
        threshold: i32,
    ) -> IndexMap<String, HashSet<String>> {
        println!{"{:?}", umis};
        // let mut adj_list: HashMap<String, HashSet<String>> = HashMap::new();
        let mut adj_list: IndexMap<String, HashSet<String>> = IndexMap::new();
        // placeholder lazy clone for the sake of completion
        // overhaul this later
        let umi_length = umis[0].len();
        let substr_idx =
            Processor::build_substr_idx(umis.clone(), umi_length.try_into().unwrap(), threshold);
        let iter_umi_pairs = Processor::iter_nearest_neighbors(umis, substr_idx)
            .unwrap()
            .neighbors;

        // logic here is screwed
        for (umi1, umi2) in &iter_umi_pairs {
            if edit_distance(&umi1, &umi2) <= threshold {
                if counts.get(umi1).unwrap() >= &((counts.get(umi2).unwrap() * 2) - 1) {
                    adj_list.entry(umi1.to_string()).and_modify(|e| {
                        e.insert(umi2.to_string());
                    });
                    adj_list.entry(umi2.to_string()).and_modify(|e| {
                        e.insert(umi2.to_string());
                    });
                }
            }
        }
        // let adj_list = Vec::from_iter(adj_list);
        println!{"{:?}", adj_list};
        return adj_list;
    }

    pub fn get_connected_components(
        &self,
        umis: Vec<String>,
        // counts: HashMap<String, i32>,
        adj_list: IndexMap<String, HashSet<String>>,
    ) -> Option<Vec<VecDeque<String>>> {
        // let mut components: IndexMap<String, HashSet<String>> = IndexMap::new();
        let mut components: Vec<VecDeque<String>> = Vec::new();
        let mut found: Vec<String> = Vec::new();

        if !adj_list.is_empty() {
            for node in adj_list.keys() {
                if !found.contains(node) {
                    let component = breadth_first_search(node.to_string(), &adj_list);
                    found.extend(component.clone());
                    components.push(component);
                }
                return Some(components);
            }
        }
        return None;
    }

    pub fn group_directional(&self, clusters: Vec<VecDeque<String>>) -> Vec<Vec<String>> {
        let mut observed: Vec<String> = Vec::new();
        let mut groups: Vec<Vec<String>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                groups.push(cluster.clone().into());
                observed.push(cluster.get(0).unwrap().to_string())
            } else {
                let mut temp_cluster: Vec<String> = Vec::new();

                for node in cluster {
                    if !observed.contains(&node) {
                        temp_cluster.push(node.clone());
                        observed.push(node);
                    }
                    groups.push(temp_cluster.clone());
                }
            }
        }
        return groups;
    }

    // see if this approach works for reverse sorting
    // https://stackoverflow.com/questions/60916194/how-to-sort-a-vector-in-descending-order-in-rust
    // adj_list.sort_by_key(|x| counts.get(&x.0));
    // adj_list.reverse();
    //
    pub fn main_grouper(&self, umis: Vec<String>, counts: HashMap<String, i32>) -> Vec<Vec<String>> {
        
        let len_umis = &umis[0].len();

        let adj_list = self.get_adj_list_directional(umis.clone(), counts, 1);
        println!{"{:?}", adj_list};
        let clusters = self.get_connected_components(umis, adj_list).unwrap();
        let final_umis = self.group_directional(clusters);

        return final_umis
    


    }
}
