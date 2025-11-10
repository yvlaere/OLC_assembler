/// Tip trimming module
/// a tip node is a node that has no incoming or no outgoing edge
/// a mergeable node is a node with exactly one incoming and one outgoing edge
/// tip trimming procedure:
/// 1. find tip nodes (indegree == 0) (ignore nodes with outdegree == 0, they will be handled because their reverse-complement is a tip)
/// 2. extend the tip node with mergeable nodes
/// 3. remove the chain if it is shorter than n nodes

use std::collections::{HashSet};
use crate::create_overlap_graph::OverlapGraph;
use crate::utils;

/// Enum for node classification
#[derive(PartialEq)]
enum NodeType {
    Mergeable = 0,   // node can be merged into a linear chain (single incoming and single outgoing edge); indegree = 1 && outdegree = 1
    Tip = 1,         // no incoming edges (graph tip); indegree = 0
    Other = 2,       // other node types
}

/// Return the list of outgoing targets for node n that currently exist in the graph.
/// (This filters out edges that point to missing nodes.)
fn target_nodes(graph: &OverlapGraph, n: &str) -> Vec<String> {
    if let Some(node) = graph.nodes.get(n) {
        node.edges.iter().map(|e| e.target_id.clone()).filter(|tgt| graph.nodes.contains_key(tgt)).collect()
    } 
    else {
        Vec::new()
    }
}

/// Classify nodes
fn node_classification(graph: &OverlapGraph, n: &str) -> (NodeType, Option<String>) {
    
    // count incoming edges to n by checking outgoing edges of rc(n)
    let incoming = target_nodes(graph, &utils::rc_node(n));
    let num_in = incoming.len();
    let outgoing = target_nodes(graph, &n);
    let num_out = outgoing.len();
    if num_in == 0 {
        return (NodeType::Tip, None);
    }
    if num_in == 1 && num_out == 1 {
        return (NodeType::Mergeable, Some(outgoing[0].clone()));
    }
    else {
        return (NodeType::Other, None);
    }
}

/// Extend a tip node up to max_ext steps
/// - collects visited nodes into chain (first entry is the tip node)
/// - returns the NodeType of the termination node (Mergeable if we reached max_ext, otherwise the non-mergeable type)
/// - returns the chain vector (the sequence of visited nodes)
fn extend(graph: &OverlapGraph, start_n: &str, max_ext: usize) -> (NodeType, Vec<String>) {
    
    // initialize
    let mut chain: Vec<String> = Vec::new();
    chain.push(start_n.to_string());
    let mut n = start_n.to_string();
    let mut steps_left = max_ext;

    // loop instead of while to guarantee a return value
    loop {

        // classify current node
        let (node_type, next_opt) = node_classification(graph, &n);

        // non-mergeable -> return
        if node_type != NodeType::Mergeable {
            return (node_type, chain);
        }

        // reached max_ext -> treat as Mergeable
        if steps_left == 0 {
            return (NodeType::Mergeable, chain);
        }

        // get next node, next_opt is an Option, so it needs to be handled
        let next = match next_opt {
            Some(s) => s,
            None => return (NodeType::Mergeable, chain),
        };

        // advance
        chain.push(next.clone());
        n = next;
        steps_left -= 1;
    }
}

/// Delete a set of nodes (both orientations) from the graph and remove associated edges
fn delete_nodes_and_edges(graph: &mut OverlapGraph, nodes_to_delete: &HashSet<String>) {
    
    // Initialize set of nodes to remove
    let mut oriented_nodes_to_delete: HashSet<String> = HashSet::new();

    for node in nodes_to_delete.iter() {
        
        // add to set of nodes to delete
        oriented_nodes_to_delete.insert(node.clone());
        // add rc counterpart
        oriented_nodes_to_delete.insert(utils::rc_node(node));
    }

    // Delete nodes from graph.nodes
    for oriented_node in oriented_nodes_to_delete.iter() {
        graph.nodes.remove(oriented_node);
    }

    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges.retain(|e| !oriented_nodes_to_delete.contains(&e.target_id));
    }
}

/// tip trimming: remove any tip nodes and their reverse-complements from the graph.
pub fn trim_tips(graph: &mut OverlapGraph, max_ext: usize) {

    // initialize collection of nodes to delete
    let mut to_delete: HashSet<String> = HashSet::new();

    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for n in keys.into_iter() {

        // skip nodes that may already be deleted
        if !graph.nodes.contains_key(&n) {
            continue;
        }

        // check if n is a Tip (only consider tips)
        let (tip_type, next) = node_classification(graph, &n);
        if tip_type != NodeType::Tip {
            continue;
        }

        // try to extend from n
        let (ext_type, chain) = extend(graph, &n, max_ext);
        // if extend returned Mergeable, skip deletion (chain may be long, not a short tip)
        if ext_type == NodeType::Mergeable {
            continue;
        }

        // otherwise the chain is small/terminating -> mark chain nodes for deletion
        for node in chain.into_iter() {
            to_delete.insert(node);
        }
    }

    if !to_delete.is_empty() {
        delete_nodes_and_edges(graph, &to_delete);
    }
}