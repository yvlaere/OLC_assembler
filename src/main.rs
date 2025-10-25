mod filter_paf;
mod create_string_graph;
mod transitive_edge_reduction;

fn main() {

    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        print_usage_and_exit(&args[0]);
    }

    let paf = &args[1];
    let max_overhang = args[2].parse::<u32>().expect("parse max_overhang");
    let overhang_ratio = args[3].parse::<f64>().expect("parse overhang_ratio");

    let (mut graph, stats) = create_string_graph(paf, max_overhang, overhang_ratio)?;
    // sort edges in each node (optional)
    //graph.sort_all_edges();

    eprintln!("Stats: internal={}, first_contained={}, second_contained={}, proper={}",
        stats.nr_internal_match, stats.nr_first_contained, stats.nr_second_contained, stats.nr_proper_overlaps);
    eprintln!("Graph nodes: {}", graph.nodes.len());

    // small demo: print degree of first N nodes
    for (i, (node_id, node)) in graph.nodes.iter().take(20).enumerate() {
        println!("{}: outdeg={} edges={:?}", node_id, node.edges.len(), node.edges.iter().take(8).collect::<Vec<_>>());
        if i >= 19 { break; }
    }

    Ok(())
}