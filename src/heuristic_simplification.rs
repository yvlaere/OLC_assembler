use crate::create_overlap_graph::OverlapGraph;

/// Remove short edges (edges with overlap length below threshold)
/// For each node with multiple outgoing edges, keep only those with overlap length >= drop_ratio * best_overlap_len
pub fn remove_short_edges(graph: &mut OverlapGraph, drop_ratio: f64) {
    let mut n_short = 0;

    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    
    for node_id in keys {
        // Get the outgoing edges for this node
        let edges_to_remove: Vec<String> = if let Some(node) = graph.nodes.get(&node_id) {
            // Skip if less than 2 outgoing edges
            if node.edges.len() < 2 {
                continue;
            }

            // Find the maximum overlap length
            let max_overlap = node.edges.iter().map(|e| e.overlap_len).max().unwrap_or(0);
            
            // Calculate threshold
            let threshold = (max_overlap as f64 * drop_ratio + 0.499) as u32;

            // Collect edges that are below the threshold
            node.edges.iter().filter(|e| e.overlap_len < threshold).map(|e| e.target_id.clone()).collect()
        } else {
            continue;
        };

        // Remove the short edges
        for target_id in edges_to_remove {
            if let Some(node) = graph.nodes.get_mut(&node_id) {
                node.remove_edge(&target_id);
                n_short += 1;
                // Remove the reverse edge as well
                if let Some(target_node) = graph.nodes.get_mut(&target_id) {
                    target_node.remove_edge(&node_id);
                }
            }
        }
    }

    eprintln!("remove short edges: removed {} short edges", n_short);
}

/// Remove low identity from nodes with multiple outgoing edges
pub fn remove_weak(graph: &mut OverlapGraph,) {
    
    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for n in keys.into_iter() {

        // check the amount of outgoing edges
        let outgoing = match graph.nodes.get(&n) {
            Some(n) => n.edges.iter().map(|e| (e.target_id.clone(), e.identity)).collect::<Vec<_>>(),
            None => continue,
        };

        // skip if less than 2 outgoing edges
        if outgoing.len() < 2 {
            continue;
        }
        
        // only keep edge with highest identity, remove others
        let mut max_identity: f64 = -1.0;
        let mut best_target: Option<String> = None;
        for (target_id, identity) in outgoing.iter() {
            if *identity > max_identity {
                max_identity = *identity;
                best_target = Some(target_id.clone());
            }
        }

        // remove all edges except the best one
        for (target_id, _) in outgoing.iter() {
            if Some(target_id) != best_target.as_ref() {
                if let Some(node) = graph.nodes.get_mut(&n) {
                    node.remove_edge(target_id);
                    // remove the reverse edge as well
                    if let Some(target_node) = graph.nodes.get_mut(target_id) {
                        target_node.remove_edge(&n);
                    }
                }
            }
        }
    }
}