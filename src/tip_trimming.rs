use std::collections::{HashSet};
use crate::create_overlap_graph::OverlapGraph;

/// Return codes for node classification
#[derive(Debug, PartialEq, Eq)]
enum NodeType {
    Mergeable = 0,   // node can be merged into a linear chain (single incoming and single outgoing edge) indegree = 1 && outdegree = 1
    Tip = 1,         // no incoming edges (graph tip); indegree = 0
    MultiOut = 2,    // multiple incoming edges (branching into this node); indegree > 1
    MultiNei = 3,    // the unique predecessor's outgoing degree != 1 (neighbor is branching); indegree = 1 && predecessor outdegree != 1
}

/// Reverse-complement (flip trailing '+' <-> '-').
fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len()-1];
            return format!("{}-", base);
        } else if last == '-' {
            let base = &id[..id.len()-1];
            return format!("{}+", base);
        }
    }
    id.to_string()
}

/// Return the list of outgoing targets for node n that currently exist in the graph.
/// (This filters out edges that point to missing nodes.)
fn out_targets(graph: &OverlapGraph, n: &str) -> Vec<String> {
    if let Some(node) = graph.nodes.get(n) {
        node.edges
            .iter()
            .map(|e| e.target_id.clone())
            .filter(|tgt| graph.nodes.contains_key(tgt))
            .collect()
    } 
    else {
        Vec::new()
    }
}

/// Degree-based classification of a node
///
/// Returns (NodeType, Option<next_node>) where next_node is valid only when returned NodeType == Mergeable.
/// The "next_node" corresponds to the single outgoing target when the local pattern is mergeable.
fn node_classification(graph: &OverlapGraph, n: &str) -> (NodeType, Option<String>) {
    
    // Count incoming edges to `n` by checking outgoing targets of rc(n)
    let incoming = out_targets(graph, &rc_node(n));
    let num_in = incoming.len();
    if num_in == 0 {
        return (NodeType::Tip, None);
    }
    if num_in > 1 {
        return (NodeType::MultiOut, None);
    }

    // indegree == 1
    // Exactly one incoming oriented node exists. `incoming[0]` is the oriented id
    let incoming_oriented = incoming.into_iter().next().unwrap();

    // now we flip the orientation to match that of n
    let w = rc_node(&incoming_oriented);

    // Outgoing targets from the predecessor-oriented node `w`.
    let outs_w = out_targets(graph, &w);
    let num_outs_w = outs_w.len();
    if num_outs_w != 1 {
        return (NodeType::MultiNei, Some(w));
    }

    // Finally, determine the single outgoing target of `n` to return as `next` when
    // mergeable. If `n` does not have exactly one outgoing target, consider it
    // non-mergeable (treat as MultiNei to be conservative).
    let outs = out_targets(graph, n);
    if outs.len() != 1 {
        return (NodeType::MultiNei, Some(w));
    }
    let next = outs.into_iter().next().unwrap();
    (NodeType::Mergeable, Some(next))
}

/// Extend along mergeable links starting from `start_n` up to `max_ext` steps.
/// - collects oriented node ids (strings) visited into `chain` (first entry is start_n),
/// - returns the termination NodeType (MERGEABLE if we used up the budget or didn't hit a non-mergeable state),
/// - and the chain vector (the sequence of oriented nodes visited).
fn extend(graph: &OverlapGraph, start_n: &str, max_ext: usize) -> (NodeType, Vec<String>) {
    let mut chain: Vec<String> = Vec::new();
    chain.push(start_n.to_string());
    let mut n = start_n.to_string();
    let mut steps_left = max_ext;

    loop {
        // classify current node
        let (ret, maybe_next) = node_classification(graph, &n);
        if ret != NodeType::Mergeable {
            return (ret, chain);
        }

        // Mergeable: get next node and continue
        let next = match maybe_next {
            Some(n) => n,
            None => return (NodeType::Mergeable, chain), // defensive; should not happen when Mergeable
        };

        // reached budget -> treat as MERGEABLE (do not delete)
        if steps_left == 0 {
            return (NodeType::Mergeable, chain);
        }

        // advance
        chain.push(next.clone());
        n = next;
        steps_left -= 1;
    }
}

/// Delete a set of sequence IDs (both orientations) from the graph and remove incident edges.
/// Input `to_delete_seqs` are oriented node ids. The whole sequence (both + and -) gets deleted.
/// Returns number of nodes (oriented ids) actually removed (counting both orientations).
fn delete_vertices_and_cleanup(graph: &mut OverlapGraph, to_delete_oriented: &HashSet<String>) -> usize {
    
    // Initialize set of oriented nodes to remove
    let mut oriented_to_remove: HashSet<String> = HashSet::new();

    for oriented in to_delete_oriented.iter() {
        
        // add to set of nodes to delete
        oriented_to_remove.insert(oriented.clone());
        // add rc counterpart
        oriented_to_remove.insert(rc_node(oriented));
    }

    // Delete nodes from graph.nodes
    for oid in oriented_to_remove.iter() {
        graph.nodes.remove(oid);
    }

    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges.retain(|e| !oriented_to_remove.contains(&e.target_id));
    }

    oriented_to_remove.len()
}

/// Delete a single directed edge src -> tgt if present.
/// Returns true if an edge was removed.
fn delete_edge(graph: &mut OverlapGraph, src: &str, tgt: &str) -> bool {
    if let Some(node) = graph.nodes.get_mut(src) {
        let before = node.edges.len();
        node.edges.retain(|e| e.target_id != tgt);
        return node.edges.len() != before;
    }
    false
}

/// tip cutting: remove any tip nodes and their reverse-complements from the graph.
/// - max_ext: extension budget in number of steps
/// returns number of oriented nodes removed (counting both orientations)
pub fn cut_tips(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed_count = 0usize;

    // We'll collect to-delete oriented nodes across the full scan and remove them in a batch
    let mut to_delete: HashSet<String> = HashSet::new();

    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for n in keys.into_iter() {

        // Skip nodes that may already be deleted
        if !graph.nodes.contains_key(&n) {
            continue;
        }
        // Check if n is a TIP (only consider tips)
        let (kind, _maybe_next) = node_classification(graph, &n);
        if kind != NodeType::Tip {
            continue;
        }

        // Attempt to extend from n
        let (ret, chain) = extend(graph, &n, max_ext);
        // If extend returned Mergeable, skip deletion (chain may be long, not a short tip)
        if ret == NodeType::Mergeable {
            continue;
        }

        // Otherwise the chain is small/terminating -> mark chain nodes for deletion
        for oid in chain.into_iter() {
            to_delete.insert(oid);
        }
    }

    if !to_delete.is_empty() {
        removed_count += delete_vertices_and_cleanup(graph, &to_delete);
    }

    // return number of oriented nodes removed (both + and - counted)
    removed_count
}

/*
/// cut internal:
/// remove short internal sequences that are bracketed by multi-neighbor conditions.
pub fn cut_internal(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed_count = 0usize;
    let mut to_delete: HashSet<String> = HashSet::new();

    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for n in keys.into_iter() {
        if !graph.nodes.contains_key(&n) { continue; }
        let (kind, _maybe_next) = node_classification(graph, &n);
        if kind != NodeType::MultiNei { continue; }
        let (ret, chain) = extend(graph, &n, max_ext);
        if ret != NodeType::MultiNei { continue; }
        for oid in chain.into_iter() { to_delete.insert(oid); }
    }

    if !to_delete.is_empty() {
        removed_count += delete_vertices_and_cleanup(graph, &to_delete);
    }

    removed_count
}

/// cut biloop:
/// This will remove one edge of a small bi-loop if the overlap lengths indicate one side is weaker.
pub fn cut_biloop(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed = 0usize;
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

    for n in keys.into_iter() {
        if !graph.nodes.contains_key(&n) { continue; }
        let (type_n, _mn) = node_classification(graph, &n);
        if type_n != NodeType::MultiNei { continue; }

        let (ret, chain) = extend(graph, &n, max_ext);
        if ret != NodeType::MultiOut { continue; }
        if chain.is_empty() { continue; }
        let last = chain.last().unwrap().clone();
        let x = rc_node(&last); // flip orientation for x

        // find the unique (non-deleted) outgoing neighbor w of n
        let outs_n = out_targets(graph, &n);
        if outs_n.is_empty() { continue; }
        let w = outs_n[0].clone();

        // inspect outgoing edges of w to find overlap lengths w -> n and w -> x
        let mut on: u32 = 0;
        let mut ox: u32 = 0;
            if let Some(node_w) = graph.nodes.get(&w) {
            for e in node_w.edges.iter() {
                if e.target_id == x {
                    ox = e.edge_len;
                }
                if e.target_id == n {
                    on = e.edge_len;
                }
            }
        }

        if on == 0 && ox == 0 {
            continue;
        }

        if on > ox {
            // delete the edge w->x and its reverse complement rc(x) -> rc(w)
            if delete_edge(graph, &w, &x) {
                // also delete the reverse complement edge: rc(x) -> rc(w)
                let rc_x = rc_node(&x);
                let rc_w = rc_node(&w);
                delete_edge(graph, &rc_x, &rc_w);
                removed += 1;
            }
        }
    }

    removed
}
*/