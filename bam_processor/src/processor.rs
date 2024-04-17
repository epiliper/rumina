extern crate bam;
extern crate rust_htslib;
use indexmap::IndexMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use strsim::hamming;
use std::process;


// pub struct NeighboringUMIs {
//     pub neighbors: HashMap<String, String>,
// }

// impl NeighboringUMIs {
//     pub fn iter(&self) -> NeighborIterator {
//         NeighborIterator {
//             nbrs: &self.neighbors,
//             index: 0,
//         }
//     }
// }

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

pub struct Processor<'b> {
    // pub groups: UMIReads,
    pub umis: &'b Vec<String>,
}
impl<'b> Processor<'b>{

    pub fn breadth_first_search(
        mut node: &'b String,
        adj_list: IndexMap<&'b String, HashSet<&'b String>>,
        // ) -> IndexMap<String, HashSet<String>>
    ) -> VecDeque<&'b String> {
        // let mut searched: IndexMap<String, HashSet<String>> = IndexMap::new();
        // let mut queue: IndexMap<String, HashSet<String>> = IndexMap::new();
        let mut searched: VecDeque<&String> = VecDeque::new();
        let mut queue: VecDeque<&String> = VecDeque::new();

        queue.push_back(node);

        while queue.len() > 0 {
            node = &queue.pop_front().unwrap();
            for next_node in &adj_list[node] {
                if !searched.contains(&next_node) {
                    queue.push_back(next_node);
                    searched.push_back(next_node);
                }
            }
        }
        return searched;
    }

    pub fn edit_distance(ua: &String, ub: &String) -> Result<usize, strsim::StrSimError> {
        return hamming(ua, ub);
    }

    pub fn get_adj_list_directional(
        &self,
        // mut umis: Vec<String>,
        counts: HashMap<String, i32>,
        threshold: usize,
    ) -> IndexMap<&'b String, HashSet<&'b String>> {
        println!{"{}", self.umis.len()};
        println!{"getting adjacency list..."};
        let mut adj_list: IndexMap<&'b String, HashSet<&'b String>> = IndexMap::new();
        let umi_length = self.umis[0].len();
        let mut i = 0;
        while i < self.umis.len() {
            let top = &self.umis[i];
            i += 1;
            let remainder = &self.umis[i..];
            for sub in remainder {
                if Processor::edit_distance(top, sub).unwrap() <= threshold {
                    if counts.get(top).unwrap() > &((counts.get(sub).unwrap() * 2 - 1)) {
                        adj_list.entry(top).and_modify(|e| {
                            e.insert(sub);
                        });
                        adj_list.entry(sub).and_modify(|e| {
                            e.insert(top);
                        });
                    }
                }
            }
        }

    println!{"{:?}", adj_list};
    return adj_list;
    }

    pub fn get_connected_components(
        &self,
        // umis: Vec<String>,
        // counts: HashMap<String, i32>,
        adj_list: IndexMap<&String, HashSet<&String>>,
    ) -> Option<Vec<VecDeque<String>>> {
        // let mut components: IndexMap<String, HashSet<String>> = IndexMap::new();
        println!{"getting connected components..."};
        let mut components: Vec<VecDeque<String>> = Vec::new();
        let mut found: Vec<&String> = Vec::new();

        if adj_list.len() > 0 {
            for node in adj_list.keys() {
                if !found.contains(node) {
                    let component = Processor::breadth_first_search(node, adj_list.clone());
                    found.extend(&component);
                    // components.push(component.iter().map(|x| *x.to_string()).collect::<VecDeque<_>>());
                    components.push(component.iter().map(|x| x.to_string()).collect::<VecDeque<_>>());
                }
            }
            return Some(components);
        } else {
            println! {"no dice"};
            return None;
        }
    }

    pub fn group_directional(&self, clusters: Vec<VecDeque<String>>) -> Vec<Vec<String>> {
        println!{"generating groups...."};
        let mut observed: Vec<String> = Vec::new();
        let mut groups: Vec<Vec<String>> = Vec::new();

        for cluster in clusters {
            if cluster.len() == 1 {
                groups.push(cluster.clone().into());
                observed.push(cluster.get(0).unwrap().to_string());
            } else {
                let mut temp_cluster: Vec<String> = Vec::new();

                for node in cluster {
                    if !observed.contains(&&node) {
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
    pub fn main_grouper(
        &self,
        umis: Vec<String>,
        counts: HashMap<String, i32>,
    ) -> Option<Vec<Vec<String>>> {
        // let len_umis = &umis[0].len();
        let final_umis: Vec<Vec<String>>;

        let adj_list = self.get_adj_list_directional(counts, 1);
        if adj_list.len() > 0 {
            let clusters = self.get_connected_components(adj_list).unwrap();
            let final_umis = self.group_directional(clusters);
            return Some(final_umis);

        } else {
            println!{"Insufficient Data"};
            return None
        }

    }
}
