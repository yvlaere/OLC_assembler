use crate::create_overlap_graph::OverlapGraph;

// this breaks the synchronization of the bigraph!
// multiple incoming edges are still possible

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
                }
            }
        }
    }
}