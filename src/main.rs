mod filter_paf;
mod create_overlap_graph;
mod transitive_edge_reduction;
mod graph_analysis;
mod compress_graph;
mod tip_trimming;

use std::env;
use std::io;

fn main() -> io::Result<()>{

    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <input.paf> <intermediary.paf> <max_overhang> <overhang_ratio>", args[0]);
        std::process::exit(1);
    }

    println!("Filtering PAF: {}", &args[1]);

    // First filter the PAF file
    filter_paf::filter_paf(&args[1], &args[2]);

    println!("Filtered PAF written to {}", &args[2]);

    let paf = &args[2];
    let max_overhang = args[3].parse::<u32>().expect("parse max_overhang");
    let overhang_ratio = args[4].parse::<f64>().expect("parse overhang_ratio");

    // Then create the overlap graph
    let (mut graph, stats) = create_overlap_graph::create_overlap_graph(paf, max_overhang, overhang_ratio)?;

    println!("Stats: internal={}, first_contained={}, second_contained={}, proper={}",
        stats.nr_internal_match, stats.nr_first_contained, stats.nr_second_contained, stats.nr_proper_overlaps);
    println!("Graph nodes: {}", graph.nodes.len());

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

    // graph analysis: weakly connected components
    let components = graph_analysis::weakly_connected_components(&graph);
    println!("Number of weakly connected components: {}", components.len());
    graph_analysis::component_sizes_sorted(&graph)
        .into_iter()
        .enumerate()
        .for_each(|(i, size)| {
            println!("Component {}: size {}", i + 1, size);
        });

    // Reduce transitive edges again
    let fuzz = 10;
    transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);

    //number of edges after reduction
    let after_edge_count_2: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
    println!("Graph edges after reduction: {}", after_edge_count_2);

    // difference
    println!("Number of removed edges: {}", after_edge_count - after_edge_count_2);

    // node to edge ratio
    let node_count = graph.nodes.len();
    let edge_count = after_edge_count_2;
    let node_to_edge_ratio = node_count as f64 / edge_count as f64;
    println!("Graph nodes: {}", node_count);
    println!("Node to edge ratio: {:.4}", node_to_edge_ratio);

    // Check that the bigraph is synchronized
    transitive_edge_reduction::check_synchronization(&graph);

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

    // Calculate compressible nodes (in_degree = out_degree = 1)
    let potentially_compressible: usize = graph.nodes.iter()
        .filter(|(id, node)| {
            let in_deg = indegree_dist.get(&1).copied().unwrap_or(0);
            let out_deg = node.edges.len();
            in_deg == 1 && out_deg == 1
        })
        .count();
    println!("\nNodes with in_degree = out_degree = 1 (potentially compressible): {}", potentially_compressible);
    
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

    let removed = tip_trimming::trim_tips(&mut graph, 4);
    println!("Removed {} tip nodes", removed);

    let tip_lengths = graph_analysis::tip_length_distribution(&graph, 1000);
    let short_tips = tip_lengths.iter().filter(|&&l| l <= 3).count();
    println!("Tip count (starts): {}, short tips (<=3): {}", tip_lengths.len(), short_tips);

    Ok(())
}