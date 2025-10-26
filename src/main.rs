mod filter_paf;
mod create_string_graph;
mod transitive_edge_reduction;
mod graph_analysis;

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

    // Then create the string graph
    let (mut graph, stats) = create_string_graph::create_string_graph(paf, max_overhang, overhang_ratio)?;

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

    Ok(())
}