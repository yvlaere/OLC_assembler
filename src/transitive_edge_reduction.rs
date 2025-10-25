use crate::create_string_graph::AssemblyGraph;

use std::collections::{HashMap, HashSet};

/// Enum for node marking (vacant, in-play, eliminated)
enum Mark {
    Vacant,
    InPlay,
    Eliminated,
}

/// create a node id for the reverse complement of the given node id
fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &name[..name.len()-1];
            format!("{}-", base)
        } else if last == '-' {
            let base = &name[..name.len()-1];
            format!("{}+", base)
        } else {
            eprintln!("Warning: unexpected node id format: {}", id);
        }
    } else {
        name.to_string()
    }
}

/// Reduce transitive edges. transitive edges are redundant edges that don't add any information to the graph.
/// Say read 1 overlaps with read 2 and read 2 overlaps with read 3 and read 1 also overlaps with read 3, then this last overlap is redundant, represented by a transitive edge.
/// Algorithm based on https://doi.org/10.1093/bioinformatics/bti1114
pub fn reduce_transitive_edges(g: &mut AssemblyGraph, fuzz: u32) -> usize {

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
    let mut nr_reduced: usize = 0;
    // We must collect a snapshot of nodes because we mutate graph during removal
    let node_keys2: Vec<String> = g.nodes.keys().cloned().collect();
    for n1 in node_keys2.iter() {
        // if node disappeared, skip
        if let Some(node) = g.nodes.get_mut(n1) {
            // collect targets to consider (we iterate over clone to avoid borrow conflicts)
            let targets: Vec<String> = node.edges.iter().map(|(t, _)| t.clone()).collect();
            for n2 in targets {
                if reduced.contains(&(n1.clone(), n2.clone())) {
                    // remove the edge n1 -> n2
                    node.remove_edge(&n2);
                    nr_reduced += 1;

                    // remove the reverse complement pair if present:
                    let n1_rc = rc_node(n1);
                    let n2_rc = rc_node(&n2);
                    if let Some(node_rc) = g.nodes.get_mut(&n2_rc) {
                        // check if node_rc has edge to n1_rc and remove if so
                        if node_rc.edges.iter().any(|(t, _)| t == &n1_rc) {
                            node_rc.remove_edge(&n1_rc);
                            nr_reduced += 1;
                        }
                    }
                }
            }
        }
    }

    nr_reduced
}