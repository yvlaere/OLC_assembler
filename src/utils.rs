/// General functions used across the project

use std::collections::{HashSet};
use crate::create_overlap_graph::OverlapGraph;

/// Get the reverse-complement of a node (flip trailing '+' <-> '-').
pub fn rc_node(id: &str) -> String {
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

/// Delete a set of nodes (both orientations) from the graph and remove associated edges
pub fn delete_nodes_and_edges(graph: &mut OverlapGraph, nodes_to_delete: &HashSet<String>) {
    
    // Initialize set of nodes to remove
    let mut oriented_nodes_to_delete: HashSet<String> = HashSet::new();

    for node in nodes_to_delete.iter() {
        
        // add to set of nodes to delete
        oriented_nodes_to_delete.insert(node.clone());
        // add rc counterpart
        oriented_nodes_to_delete.insert(rc_node(node));
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