mod alignment_filtering;
mod filter_paf;
mod create_overlap_graph;
mod transitive_edge_reduction;
mod graph_analysis;
mod compress_graph;
mod tip_trimming;
mod bubble_removal;
mod utils;

use std::env;
use std::io;
use std::collections::HashSet;

fn main() -> io::Result<()>{

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.paf> <intermediary.paf>", args[0]);
        std::process::exit(1);
    }

    println!("Filtering PAF: {}", &args[1]);

    // First filter the PAF file

    // Configuration
    let min_overlap_length: u32 = 2000;
    let min_overlap_count: u32 = 3;
    let min_percent_identity: f32 = 5.0; // 100/2000 = 5%
    let max_overhang = 0;
    let overhang_ratio = 0.8;
    let overlaps = alignment_filtering::filter_paf(&args[1], &args[2], &min_overlap_length, &min_overlap_count, &min_percent_identity, &max_overhang, &overhang_ratio);

    println!("Filtered PAF written to {}\n", &args[2]);

    // Then create the overlap graph
    let mut graph = create_overlap_graph::create_overlap_graph(overlaps?)?;

    // Check that the bigraph is synchronized
    graph_analysis::check_synchronization(&graph);

    // Iterative graph cleanup until convergence
    println!("=== STARTING GRAPH CLEANUP ===");
    let max_bubble_len = 8;  // Max edges per bubble path
    let min_support_ratio = 1.1;  // Require 10% stronger path for removal
    let max_tip_len = 4;  // Maximum tip length to remove
    let fuzz = 10;  // Fuzz factor for transitive reduction

    let mut iteration = 1;
    let mut prev_node_count = graph.nodes.len() + 1;  // Ensure first iteration runs

    for i in 0..10 {
        prev_node_count = graph.nodes.len();
        println!("\n=== Cleanup Iteration {} ===", iteration);

        // 1) Transitive edge reduction
        let edges_before: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
        transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);
        let edges_after: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
        println!("Removed {} edges with transitive edge reduction", edges_before.saturating_sub(edges_after));

        // 2) Bubble popping
        let node_count_before_bubble_popping = graph.nodes.len();
        bubble_removal::remove_bubbles(&mut graph, max_bubble_len, min_support_ratio);
        let node_count_after_bubble_popping = graph.nodes.len();
        println!("Removed {} bubble nodes (including RCs)", node_count_before_bubble_popping - node_count_after_bubble_popping);

        // 3) Remove small components (size < 2)
        let components = graph_analysis::weakly_connected_components(&graph);
        let mut comp_nodes_to_remove: HashSet<String> = HashSet::new();
        for component in components.iter() {
            if component.len() < 2 {
                for nid in component.iter() {
                    comp_nodes_to_remove.insert(nid.clone());
                }
            }
        }
        let small_comp_count = comp_nodes_to_remove.len();
        for node_id in comp_nodes_to_remove.iter() {
            graph.nodes.remove(node_id.as_str());
            // also remove RC oriented node if present
            if node_id.ends_with('+') {
                let rc = node_id[..node_id.len()-1].to_string() + "-";
                graph.nodes.remove(rc.as_str());
            } else if node_id.ends_with('-') {
                let rc = node_id[..node_id.len()-1].to_string() + "+";
                graph.nodes.remove(rc.as_str());
            }
        }
        println!("Removed {} oriented nodes from small components (<2)", small_comp_count);

        // 4) Tip trimming
        let node_count_before_trimming = graph.nodes.len();
        tip_trimming::trim_tips(&mut graph, max_tip_len);
        let node_count_after_trimming = graph.nodes.len();
        println!("Removed {} nodes by tip trimming", node_count_before_trimming - node_count_after_trimming);

        // 5) Iteration stats
        let (_indegree_dist, outdegree_dist) = graph_analysis::analyze_degrees(&graph);
        let high_branch = outdegree_dist.iter().filter(|(deg, _)| **deg > 1).map(|(_, count)| *count).sum::<usize>();
        let very_high = outdegree_dist.iter().filter(|(deg, _)| **deg >= 3).map(|(_, count)| *count).sum::<usize>();

        println!("\nIteration {} stats:", iteration);
        println!("Current nodes: {}", graph.nodes.len());
        let current_edges: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
        println!("Current edges: {}", current_edges);
        println!("Node/edge ratio: {:.4}", graph.nodes.len() as f64 / current_edges as f64);
        println!("Nodes with out-degree > 1: {}", high_branch);
        println!("Nodes with out-degree >= 3: {}", very_high);

        iteration += 1;

        println!("\n");
        println!("\n----------------------------------------");
        println!("\n");
    }

    // weakly connected components after cleanup
    let final_components = graph_analysis::weakly_connected_components(&graph);
    println!("Final weakly connected components: {}", final_components.len());

    // graph compression into unitigs
    println!("Compressing graph into unitigs...");
    let compressed_graph = compress_graph::compress_unitigs(&graph);
    println!("Compressed graph has {} unitigs", compressed_graph.unitigs.len());
    
    // unitig length stats
    let mut unitig_lengths: Vec<usize> = compressed_graph.unitigs.iter()
        .map(|u| u.members.len())
        .collect();
    unitig_lengths.sort_unstable();
    let total_unitigs = unitig_lengths.len();
    let total_length: usize = unitig_lengths.iter().sum();
    let N50 = {
        let mut cum_length = 0;
        let half_length = total_length / 2;
        let mut n50 = 0;
        for &len in unitig_lengths.iter().rev() {
            cum_length += len;
            if cum_length >= half_length {
                n50 = len;
                break;
            }
        }
        n50
    };
    println!("Unitig length stats:");
    println!("Total unitigs: {}", total_unitigs);
    println!("Total length (in nodes): {}", total_length);
    println!("N50 unitig length (in nodes): {}", N50);

    // largest unitig members
    if let Some(largest_unitig) = compressed_graph.unitigs.iter().max_by_key(|u| u.members.len()) {
        println!("Largest unitig ID: {}, length (in nodes): {}", largest_unitig.id, largest_unitig.members.len());
    }

    // Final stats
    println!("\n=== Final Graph Stats ===");
    let final_node_count = graph.nodes.len();
    let final_edge_count: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
    println!("Final graph nodes: {}", final_node_count);
    println!("Final graph edges: {}", final_edge_count);
    println!("Final node to edge ratio: {:.4}", final_node_count as f64 / final_edge_count as f64);

    Ok(())
}