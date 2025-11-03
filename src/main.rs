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
    println!("Graph nodes: {}\n", graph.nodes.len());

    //number of edges before reduction
    let initial_edge_count: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
    println!("Graph edges before reduction: {}", initial_edge_count);

    // node to edge ratio
    let node_to_edge_ratio = graph.nodes.len() as f64 / initial_edge_count as f64;
    println!("Node to edge ratio: {:.4}", node_to_edge_ratio);

    // Check that the bigraph is synchronized
    transitive_edge_reduction::check_synchronization(&graph);

    // Reduce transitive edges
    let fuzz = 10;
    transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);

    //number of edges after reduction
    let after_edge_count: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
    println!("Graph edges after reduction: {}", after_edge_count);

    // difference
    println!("Number of removed edges: {}", initial_edge_count - after_edge_count);

    // node to edge ratio
    let node_count = graph.nodes.len();
    let edge_count = after_edge_count;
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
        let removed_tips = tip_trimming::trim_tips(&mut graph, max_tip_len);
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
    }

    println!("\nGraph cleanup converged after {} iterations", iteration - 1);

    let tip_lengths = graph_analysis::tip_length_distribution(&graph, 1000);
    let short_tips = tip_lengths.iter().filter(|&&l| l <= 3).count();
    println!("Tip count (starts): {}, short tips (<=3): {}", tip_lengths.len(), short_tips);

    // Analyze graph before bubble removal
    println!("\nAnalyzing graph before bubble removal:");
    let initial_nodes = graph.nodes.len();
    let (indegree_dist, outdegree_dist) = graph_analysis::analyze_degrees(&graph);
    
    // Count branching nodes (in or out > 1)
    let branching = outdegree_dist.iter()
        .filter(|(deg, _)| **deg > 1)
        .map(|(_, count)| *count)
        .sum::<usize>();
    println!("Nodes with out-degree > 1: {}", branching);
    
    // High branching (>= 3 out edges)
    let high_branch = outdegree_dist.iter()
        .filter(|(deg, _)| **deg >= 3)
        .map(|(_, count)| *count)
        .sum::<usize>();
    println!("Nodes with out-degree >= 3: {}", high_branch);

    // Bubble popping already performed in the iterative cleanup above.

    // Re-analyze graph structure
    let (indegree_dist, outdegree_dist) = graph_analysis::analyze_degrees(&graph);
    
    // Count remaining branching
    let remaining_branch = outdegree_dist.iter()
        .filter(|(deg, _)| **deg > 1)
        .map(|(_, count)| *count)
        .sum::<usize>();
    println!("Remaining nodes with out-degree > 1: {}", remaining_branch);
    
    let remaining_high = outdegree_dist.iter()
        .filter(|(deg, _)| **deg >= 3)
        .map(|(_, count)| *count)
        .sum::<usize>();
    println!("Remaining nodes with out-degree >= 3: {}", remaining_high);
    
    // Check synchronization after bubble removal
    transitive_edge_reduction::check_synchronization(&graph);

    // node to edge ratio
    let node_count = graph.nodes.len();
    let edge_count = after_edge_count;
    let node_to_edge_ratio = node_count as f64 / edge_count as f64;
    println!("Graph nodes: {}", node_count);
    println!("Node to edge ratio: {:.4}", node_to_edge_ratio);

    // graph analysis: weakly connected components
    let components = graph_analysis::weakly_connected_components(&graph);
    println!("Number of weakly connected components: {}", components.len());
    graph_analysis::component_sizes_sorted(&graph)
        .into_iter()
        .enumerate()
        .for_each(|(i, size)| {
            println!("Component {}: size {}", i + 1, size);
        });

    // Analyze graph structure before compression
    println!("\nAnalyzing graph structure before compression:");
    let (indegree_dist, outdegree_dist) = graph_analysis::analyze_degrees(&graph);
    
    // Print degree distributions
    println!("\nIn-degree distribution (degree -> count):");
    let mut degrees: Vec<_> = indegree_dist.keys().collect();
    degrees.sort_unstable();
    for &deg in &degrees {
        println!("{}: {}", deg, indegree_dist[&deg]);
    }
    
    println!("\nOut-degree distribution (degree -> count):");
    let mut degrees: Vec<_> = outdegree_dist.keys().collect();
    degrees.sort_unstable();
    for &deg in &degrees {
        println!("{}: {}", deg, outdegree_dist[&deg]);
    }

    // Note: compute-oriented compressible stats via graph_analysis::compressible_node_stats()
    
    // Create compressed graph
    let compressed_graph = compress_graph::compress_unitigs(&graph);
    println!("Number of unitigs: {}", compressed_graph.unitigs.len());
    
    // Report compression statistics
    let compression_ratio = (graph.nodes.len() - compressed_graph.unitigs.len()) as f64 / graph.nodes.len() as f64;
    println!("Compression ratio: {:.2}%", compression_ratio * 100.0);
    
    // Report average unitig length
    let avg_unitig_len = compressed_graph.unitigs.iter()
        .map(|u| u.members.len())
        .sum::<usize>() as f64 / compressed_graph.unitigs.len() as f64;
    println!("Average unitig length: {:.2}", avg_unitig_len);

    let (compressible, total, frac) = graph_analysis::compressible_node_stats(&graph);
    println!("Compressible oriented nodes: {} / {} ({:.2}%)", compressible, total, frac * 100.0);

    let tip_lengths = graph_analysis::tip_length_distribution(&graph, 1000);
    let short_tips = tip_lengths.iter().filter(|&&l| l <= 3).count();
    println!("Tip count (starts): {}, short tips (<=3): {}", tip_lengths.len(), short_tips);

    let top_branch = graph_analysis::branching_summary(&graph, 20);
    println!("Top branching nodes (id, in, out):");
    for (id, in_d, out_d) in top_branch {
        println!("  {}: in={}, out={}", id, in_d, out_d);
    }

    

    Ok(())
}