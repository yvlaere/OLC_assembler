mod filter_paf;
mod create_string_graph;
mod transitive_edge_reduction;

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

    eprintln!("Stats: internal={}, first_contained={}, second_contained={}, proper={}",
        stats.nr_internal_match, stats.nr_first_contained, stats.nr_second_contained, stats.nr_proper_overlaps);
    eprintln!("Graph nodes: {}", graph.nodes.len());

    // small demo: print degree of first N nodes
    for (i, (node_id, node)) in graph.nodes.iter().take(20).enumerate() {
        println!("{}: outdeg={} edges={:?}", node_id, node.edges.len(), node.edges.iter().take(8).collect::<Vec<_>>());
        if i >= 19 { break; }
    }

    // Finally reduce transitive edges
    let fuzz = 10;
    let nr_reduced = transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);
    println!("Reduced {} transitive edges", nr_reduced);

    Ok(())
}