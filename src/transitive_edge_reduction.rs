use crate::create_overlap_graph::OverlapGraph;

use std::collections::{HashMap, HashSet};

/// Enum for node marking (vacant, in-play, eliminated)
#[derive(Copy, Clone, PartialEq)]
enum Mark {
    Vacant,
    InPlay,
    Eliminated,
}

/// create a node id for the reverse complement of the given node id
fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len()-1];
            format!("{}-", base)
        } else if last == '-' {
            let base = &id[..id.len()-1];
            format!("{}+", base)
        } else {
            eprintln!("Warning: unexpected node id format: {}", id);
            id.to_string()
        }
    } else {
        id.to_string()
    }
}

/// Check if the bigraph is synchronized:
/// 1. Every node has a reverse complement.
/// 2. Ingoing edges of every node correspond to outgoing edges of its reverse complement.
pub fn check_synchronization(g: &OverlapGraph) {
    for n in g.nodes.keys() {
        // compute reverse complement node
        let n_rc = rc_node(n);

        // check that the reverse complement exists
        assert!(
            g.nodes.contains_key(&n_rc),
            "Reverse complement not found for {}. The bigraph is not synchronized.",
            n
        );

        // check that every outgoing edge has a counterpart in the reverse complement node
        if let Some(node) = g.nodes.get(n) {
            for (t, _) in &node.edges {
                let t_rc = rc_node(t);

                // get the reverse complement node for t_rc
                if let Some(rc_node) = g.nodes.get(&t_rc) {
                    // Check if there's a matching edge from t_rc to n_rc
                    assert!(
                        rc_node.edges.iter().any(|(target, _)| target == &n_rc),
                        "Corresponding edge not found for {}. The bigraph is not synchronized.",
                        n
                    );
                } else {
                    panic!(
                        "Reverse complement node {} missing for target {}. The bigraph is not synchronized.",
                        t_rc, t
                    );
                }
            }
        }
    }
}

/// Reduce transitive edges. transitive edges are redundant edges that don't add any information to the graph.
/// Say read 1 overlaps with read 2 and read 2 overlaps with read 3 and read 1 also overlaps with read 3, then this last overlap is redundant, represented by a transitive edge.
/// Algorithm based on https://doi.org/10.1093/bioinformatics/bti1114
pub fn reduce_transitive_edges(g: &mut OverlapGraph, fuzz: u32) {

    // Prepare node list to iterate deterministically and avoid borrow conflicts
    let node_keys: Vec<String> = g.nodes.keys().cloned().collect();

    // mark: per-node status (Vacant/InPlay/Eliminated)
    let mut mark: HashMap<String, Mark> = HashMap::with_capacity(g.nodes.len());
    
    // reduced set: node pairs (from, to) that should be removed
    let mut reduced: HashSet<(String, String)> = HashSet::new();

    // initialize all nodes as vacant and ensure edges are sorted ascending by length
    for n in &node_keys {
        mark.insert(n.clone(), Mark::Vacant);
        if let Some(node) = g.nodes.get_mut(n) {
            node.sort_edges();
        }
    }

    // main loop: For every node compare nodes encountered two steps into the future with those encountered one step into the future
    for n1 in &node_keys {

        // skip if node not present (may have been removed) or no outgoing edges
        let out_edges = match g.nodes.get(n1) {
            Some(node) => &node.edges,
            None => continue,
        };
        if out_edges.is_empty() {
            continue;
        }

        // 1) mark all direct neighbors of n1 as InPlay
        for (n2, _) in out_edges.iter() {
            mark.insert(n2.clone(), Mark::InPlay);
        }

        // 2) compute longest outgoing edge length from n1 (last after sort) + fuzz
        let longest = {
            let last_len = out_edges.last().unwrap().1 as u64;
            last_len + fuzz as u64
        };

        // 3) For each n2 (outgoing from n1), check n2->n3 edges
        for (n2, len_n1n2) in out_edges.iter() {

            // skip n2 if not InPlay
            if mark.get(n2).copied() != Some(Mark::InPlay) {
                continue;
            }

            // get node n2
            if let Some(node2) = g.nodes.get(n2) {
                for (n3, len_n2n3) in node2.edges.iter() {
                    // if path length n1->n2->n3 <= longest then candidate for elimination
                    let path_len = *len_n2n3 as u64 + *len_n1n2 as u64;
                    if path_len <= longest {
                        if mark.get(n3).copied() == Some(Mark::InPlay) {
                            mark.insert(n3.clone(), Mark::Eliminated);
                        }
                    }
                }
            }
        }

        // 4) Additional rule: if n2->n3 is very small (< fuzz) or is the smallest outgoing edge of n2,
        // then eliminate n3 if it is InPlay.
        for (n2, _) in out_edges.iter() {
            if let Some(node2) = g.nodes.get(n2) {
                // find min outgoing length for n2, if any
                let min_len_opt = node2.edges.iter().map(|(_, l)| *l).min();

                for (n3, len_n2n3) in node2.edges.iter() {
                    let do_eliminate = if *len_n2n3 < fuzz {
                        true
                    } else if let Some(min_len) = min_len_opt {
                        *len_n2n3 == min_len
                    } else {
                        false
                    };

                    if do_eliminate && mark.get(n3).copied() == Some(Mark::InPlay) {
                        mark.insert(n3.clone(), Mark::Eliminated);
                    }
                }
            }
        }

        // 5) Mark edges from n1 to eliminated nodes for removal, then reset marks of the direct neighbors to vacant
        for (n2, _) in out_edges.iter() {
            if mark.get(n2).copied() == Some(Mark::Eliminated) {
                reduced.insert((n1.clone(), n2.clone()));
            }
            // reset mark back to vacant for next iteration
            mark.insert(n2.clone(), Mark::Vacant);
        }
    } // end for n1

    // 6) Actually remove reduced edges from the graph (and remove reverse-complement counterparts if present)
    
    // Collect all edges to remove first
    let mut edges_to_remove: Vec<(String, String)> = Vec::new();
    
    for n1 in g.nodes.keys() {
        if let Some(node) = g.nodes.get(n1) {
            for (n2, _) in &node.edges {
                if reduced.contains(&(n1.clone(), n2.clone())) {
                    // Add both the edge and its reverse complement
                    edges_to_remove.push((n1.clone(), n2.clone()));
                    edges_to_remove.push((rc_node(n2), rc_node(n1)));
                }
            }
        }
    }

    // Now remove all edges in a separate pass
    for (from, to) in edges_to_remove {
        if let Some(node) = g.nodes.get_mut(&from) {
            node.remove_edge(&to);
        }
    }
}