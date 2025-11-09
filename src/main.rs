mod filter_paf;
mod create_overlap_graph;
mod transitive_edge_reduction;
mod graph_analysis;
mod compress_graph;
mod tip_trimming;
mod bubble_removal;

use std::env;
use std::io;
use std::collections::HashSet;

fn main() -> io::Result<()>{

    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <input.paf> <intermediary.paf> <max_overhang> <overhang_ratio>", args[0]);
        std::process::exit(1);
    }

    println!("Filtering PAF: {}", &args[1]);

    // First filter the PAF file (ignore Result; errors will surface later)
    //let _ = filter_paf::filter_paf(&args[1], &args[2]);

    println!("Filtered PAF written to {}\n", &args[2]);

    let paf = &args[2];
    let max_overhang = args[3].parse::<u32>().expect("parse max_overhang");
    let overhang_ratio = args[4].parse::<f64>().expect("parse overhang_ratio");

    // Then create the overlap graph
    let (mut graph, stats) = create_overlap_graph::create_overlap_graph(paf, max_overhang, overhang_ratio)?;

    println!("Stats: internal={}, first_contained={}, second_contained={}, proper={}",
        stats.nr_internal_match, stats.nr_first_contained, stats.nr_second_contained, stats.nr_proper_overlaps);

    //number of edges before reduction
    let edge_count: usize = graph.nodes.values().map(|n| n.edges.len()).sum();

    // node to edge ratio
    let node_count = graph.nodes.len();
    let node_to_edge_ratio = node_count as f64 / edge_count as f64;
    println!("Graph nodes: {}", node_count);
    println!("Node to edge ratio: {:.4}", node_to_edge_ratio);

    // Check that the bigraph is synchronized
    transitive_edge_reduction::check_synchronization(&graph);

    // Iterative graph cleanup until convergence
    println!("\nStarting iterative graph cleanup...");
    let max_bubble_len = 8;  // Max edges per bubble path
    let min_support_ratio = 1.1;  // Require 10% stronger path for removal
    let max_tip_len = 4;  // Maximum tip length to remove
    let fuzz = 10;  // Fuzz factor for transitive reduction

    let mut iteration = 1;
    let mut prev_node_count = graph.nodes.len() + 1;  // Ensure first iteration runs

    while graph.nodes.len() < prev_node_count {
        prev_node_count = graph.nodes.len();
        println!("\n=== Cleanup Iteration {} ===", iteration);

        // 1) Transitive edge reduction
        let edges_before: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
        transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);
        let edges_after: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
        println!("Transitive reduction removed {} edges", edges_before.saturating_sub(edges_after));

        // 2) Bubble popping
        let removed_bubbles = bubble_removal::remove_bubbles(&mut graph, max_bubble_len, min_support_ratio);
        println!("Removed {} bubble nodes (including RCs)", removed_bubbles);

        // 3) Remove small components (size < 4)
        let components = graph_analysis::weakly_connected_components(&graph);
        let mut comp_nodes_to_remove: HashSet<String> = HashSet::new();
        for component in components.iter() {
            if component.len() < 4 {
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
        println!("Removed {} oriented nodes from small components (<4)", small_comp_count);

        // 4) Tip trimming
        let tip_lengths = graph_analysis::tip_length_distribution(&graph, 1000);
        let short_tips = tip_lengths.iter().filter(|&&l| l <= 3).count();
        println!("Tips before trimming: {}, short tips (<=3): {}", tip_lengths.len(), short_tips);
        let removed_tips = tip_trimming::cut_tips(&mut graph, max_tip_len);
        println!("Removed {} tip nodes", removed_tips);

        // 5) Iteration stats
        let (_indegree_dist, outdegree_dist) = graph_analysis::analyze_degrees(&graph);
        let high_branch = outdegree_dist.iter()
            .filter(|(deg, _)| **deg > 1)
            .map(|(_, count)| *count)
            .sum::<usize>();
        let very_high = outdegree_dist.iter()
            .filter(|(deg, _)| **deg >= 3)
            .map(|(_, count)| *count)
            .sum::<usize>();

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
    

    Ok(())
}