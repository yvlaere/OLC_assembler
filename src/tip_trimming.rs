use std::collections::{HashMap, HashSet};
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
            .map(|(tgt, _len)| tgt.clone())
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

/// Extend along mergeable links starting from `start_v` up to `max_ext` steps.
/// Mirrors miniasm's asg_extend:
/// - collects oriented node ids (strings) visited into `chain` (first entry is start_v),
/// - returns the termination NodeType (MERGEABLE if we used up the budget or didn't hit a non-mergeable state),
/// - and the chain vector (the sequence of oriented vertices visited).
fn asg_extend(graph: &OverlapGraph, start_v: &str, max_ext: usize) -> (NodeType, Vec<String>) {
    let mut chain: Vec<String> = Vec::new();
    chain.push(start_v.to_string());
    let mut v = start_v.to_string();
    let mut steps_left = max_ext;

    loop {
        let (ret, maybe_next) = node_classification(graph, &v);
        if ret != NodeType::Mergeable {
            return (ret, chain);
        }
        // Mergeable: get next node and continue
        let next = match maybe_next {
            Some(n) => n,
            None => return (NodeType::Mergeable, chain), // defensive; should not happen when Mergeable
        };
        if steps_left == 0 {
            // reached budget -> treat as MERGEABLE (do not delete)
            return (NodeType::Mergeable, chain);
        }
        chain.push(next.clone());
        v = next;
        steps_left -= 1;
    }
}

/// Delete a set of sequence IDs (both orientations) from the graph and remove incident edges.
/// Input `to_delete_seqs` are *oriented* node ids; we delete the whole sequence (both + and -).
/// Returns number of nodes (oriented ids) actually removed (counting both orientations).
fn delete_vertices_and_cleanup(graph: &mut OverlapGraph, to_delete_oriented: &HashSet<String>) -> usize {
    // Expand to include rc counterparts and compute sequence base ids to delete
    let mut seq_bases_to_delete: HashSet<String> = HashSet::new();

    for oriented in to_delete_oriented.iter() {
        // strip trailing orientation char to get base (if present)
        if oriented.ends_with('+') || oriented.ends_with('-') {
            let base = oriented[..oriented.len()-1].to_string();
            seq_bases_to_delete.insert(base);
        } else {
            // fallback: treat full id as base
            seq_bases_to_delete.insert(oriented.clone());
        }
    }

    // Build set of oriented ids to remove (both + and - for each base)
    let mut oriented_to_remove: HashSet<String> = HashSet::new();
    for base in seq_bases_to_delete.iter() {
        oriented_to_remove.insert(format!("{}+", base));
        oriented_to_remove.insert(format!("{}-", base));
    }

    // Delete nodes from graph.nodes
    for oid in oriented_to_remove.iter() {
        graph.nodes.remove(oid);
    }

    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges.retain(|(tgt, _len)| !oriented_to_remove.contains(tgt));
    }

    oriented_to_remove.len()
}

/// Delete a single directed arc src -> tgt if present.
/// Returns true if an arc was removed.
fn delete_arc(graph: &mut OverlapGraph, src: &str, tgt: &str) -> bool {
    if let Some(node) = graph.nodes.get_mut(src) {
        let before = node.edges.len();
        node.edges.retain(|(t, _len)| t != tgt);
        return node.edges.len() != before;
    }
    false
}

/// CUT TIP: replicate miniasm's asg_cut_tip behavior.
/// - max_ext: extension budget in number of steps (use small numbers like 4)
/// returns number of oriented nodes removed (counting both orientations)
pub fn cut_tip(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed_count = 0usize;
    // We'll collect to-delete oriented nodes across the full scan and remove them in a batch (like miniasm)
    let mut to_delete: HashSet<String> = HashSet::new();

    // iterate over a snapshot of current node keys (we must not mutate while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for v in keys.into_iter() {
        // Skip nodes that may already be deleted
        if !graph.nodes.contains_key(&v) {
            continue;
        }
        // Check if v is a TIP (only consider tips)
        let (kind, _maybe_next) = node_classification(graph, &v);
        if kind != NodeType::Tip {
            continue;
        }
        // Attempt to extend from v
        let (ret, chain) = asg_extend(graph, &v, max_ext);
        // If extend returned Mergeable, skip deletion (chain may be long)
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

/// CUT INTERNAL: replicate asg_cut_internal
/// remove short internal sequences that are bracketed by multi-neighbor conditions.
pub fn cut_internal(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed_count = 0usize;
    let mut to_delete: HashSet<String> = HashSet::new();

    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for v in keys.into_iter() {
        if !graph.nodes.contains_key(&v) { continue; }
        let (kind, _maybe_next) = node_classification(graph, &v);
        if kind != NodeType::MultiNei { continue; }
        let (ret, chain) = asg_extend(graph, &v, max_ext);
        if ret != NodeType::MultiNei { continue; }
        for oid in chain.into_iter() { to_delete.insert(oid); }
    }

    if !to_delete.is_empty() {
        removed_count += delete_vertices_and_cleanup(graph, &to_delete);
    }

    removed_count
}

/// CUT BILOOP: replicate asg_cut_biloop.
/// This will remove one arc of a small bi-loop if the overlap lengths indicate one side is weaker.
pub fn cut_biloop(graph: &mut OverlapGraph, max_ext: usize) -> usize {
    let mut removed = 0usize;
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

    for v in keys.into_iter() {
        if !graph.nodes.contains_key(&v) { continue; }
        let (kind_v, _mn) = node_classification(graph, &v);
        if kind_v != NodeType::MultiNei { continue; }

        let (ret, chain) = asg_extend(graph, &v, max_ext);
        if ret != NodeType::MultiOut { continue; }
        // x = last entry in chain ^ 1 (flip orientation)
        if chain.is_empty() { continue; }
        let last = chain.last().unwrap().clone();
        let x = rc_node(&last); // flip orientation for x

        // find the unique (non-deleted) outgoing neighbor w of v
        let outs_v = out_targets(graph, &v);
        if outs_v.is_empty() { continue; }
        // In miniasm the loop assumes there is exactly one outgoing non-deleted arc for v in this context.
        // We'll pick the first valid one (consistent with the original).
        let w = outs_v[0].clone();

        // inspect outgoing arcs of w to find overlap lengths w->v and w->x
        let mut ov: u32 = 0;
        let mut ox: u32 = 0;
        if let Some(node_w) = graph.nodes.get(&w) {
            for (tgt, len) in node_w.edges.iter() {
                if tgt == &x {
                    ox = *len;
                }
                if tgt == &v {
                    ov = *len;
                }
            }
        }

        if ov == 0 && ox == 0 {
            continue;
        }

        if ov > ox {
            // delete the arc w->x and its reverse complement x^1 -> w^1 (here rc_node flips orientation)
            if delete_arc(graph, &w, &x) {
                // also delete the reverse complement arc: rc(x) -> rc(w)
                let rc_x = rc_node(&x);
                let rc_w = rc_node(&w);
                delete_arc(graph, &rc_x, &rc_w);
                removed += 1;
            }
        }
    }

    if removed > 0 {
        // After arc deletions we should also purge any isolated sequences if desired.
        // For parity with miniasm, you could call a cleanup routine here that removes orphan sequences.
        // We simply return the number of removed bi-loop arcs.
    }

    removed
}