mod cli;
mod configs;
mod alignment_filtering;
mod create_overlap_graph;
mod transitive_edge_reduction;
mod tip_trimming;
mod bubble_removal;
mod graph_analysis;
mod compress_graph;
mod utils;
mod heuristic_simplification;

use std::io;
use clap::Parser;
use cli::{Cli, Commands};
use std::collections::HashSet;

// used for deserializing overlaps
use std::io::BufReader;
use std::fs::File;
use crate::alignment_filtering::Overlap;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::AlignmentFiltering(args) => {
            let config: crate::configs::AlignmentFilteringConfig = args.into();
            // run filtering and serialize overlaps to the configured output
            let out = alignment_filtering::run_alignment_filtering(
                &config.input_paf,
                &config.min_overlap_length,
                &config.min_overlap_count,
                &config.min_percent_identity,
                &config.overhang_ratio,
            )?;
            out.serialize_overlaps(&config.output_overlaps)?;
            println!("Wrote overlaps to {}", config.output_overlaps);
        }
        Commands::CreateOverlapGraph(args) => {
            let config: crate::configs::CreateOverlapGraphConfig = args.into();
            // load overlaps from binary
            let file = File::open(&config.overlaps)?;
            let reader = BufReader::new(file);
            let overlaps: HashMap<(usize, usize), Overlap> = bincode::deserialize_from(reader)?;
            let graph = create_overlap_graph::run_create_overlap_graph(overlaps)?;
            graph.write_dot("overlap_before_heuristic.dot")?;
            println!("Wrote graph to overlap_before_heuristic.dot");
        }
        Commands::CreateUnitigs(args) => {
            let config: crate::configs::UnitigConfig = args.into();
            let file = File::open(&config.overlap_graph_binary)?;
            let reader = BufReader::new(file);
            let overlaps: HashMap<(usize, usize), Overlap> = bincode::deserialize_from(reader)?;
            let graph = create_overlap_graph::run_create_overlap_graph(overlaps)?;

            let out_path = std::path::Path::new(&config.output_dir).join(format!("{}.fa", config.output_prefix));
            let out_str = out_path.to_str().ok_or("invalid output path")?;
            let compressed = compress_graph::compress_unitigs(&graph, &config.reads_fq, out_str);
            println!("Compressed graph has {} unitigs, wrote to {}", compressed.unitigs.len(), out_str);
            let gfa_path = std::path::Path::new(&config.output_dir).join(format!("{}.gfa", config.output_prefix));
            let gfa_str = gfa_path.to_str().ok_or("invalid output path")?;
            compressed.write_gfa(gfa_str)?;
            println!("Wrote GFA to {}", gfa_str);
        }
        Commands::Assemble(args) => {
            let config: crate::configs::AssembleConfig = args.into();

            // ensure output directory exists
            let out_dir = std::path::Path::new(&config.output_dir);
            std::fs::create_dir_all(out_dir)?;

            // Determine the path to overlaps: either use provided overlaps or run alignment filtering
            let overlaps_path_str = if let Some(ref overlaps_file) = config.overlaps {
                // Use provided overlaps
                println!("Using provided overlaps from {}", overlaps_file);
                overlaps_file.clone()
            } else {
                // Run alignment filtering
                let input_paf = config.input_paf.as_ref().ok_or("Either --input-paf or --overlaps must be provided")?;
                
                // write overlaps into the output directory using the chosen prefix
                let overlaps_path = out_dir.join(format!("{}.overlaps.bin", config.output_prefix));
                let overlaps_path_str = overlaps_path.to_str().ok_or("invalid output path")?.to_string();

                let out = alignment_filtering::run_alignment_filtering(
                    input_paf,
                    &config.min_overlap_length,
                    &config.min_overlap_count,
                    &config.min_percent_identity,
                    &config.overhang_ratio,
                )?;
                out.serialize_overlaps(&overlaps_path_str)?;
                println!("Wrote overlaps to {}", overlaps_path_str);
                overlaps_path_str
            };

            // load overlaps, build graph
            let file = File::open(&overlaps_path_str)?;
            let reader = BufReader::new(file);
            let overlaps: HashMap<(usize, usize), Overlap> = bincode::deserialize_from(reader)?;
            let mut graph = create_overlap_graph::run_create_overlap_graph(overlaps)?;

            // Graph simplification: iterative cleanup
            graph_analysis::check_synchronization(&graph);
            println!("=== STARTING GRAPH CLEANUP ===");
            let max_bubble_len = 8usize;
            let min_support_ratio = 1.1f64;
            let max_tip_len = 4usize;
            let fuzz = 10u32;

            for iteration in 1..=10 {
                println!("\n=== Cleanup Iteration {} ===", iteration);

                // transitive edge reduction
                let edges_before: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
                transitive_edge_reduction::reduce_transitive_edges(&mut graph, fuzz);
                let edges_after: usize = graph.nodes.values().map(|n| n.edges.len()).sum();
                println!("Removed {} edges with transitive edge reduction", edges_before.saturating_sub(edges_after));

                // bubble removal
                let node_count_before = graph.nodes.len();
                bubble_removal::remove_bubbles(&mut graph, max_bubble_len, min_support_ratio);
                let node_count_after = graph.nodes.len();
                println!("Removed {} bubble nodes (including RCs)", node_count_before.saturating_sub(node_count_after));

                // remove small components (<2)
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
                    if node_id.ends_with('+') {
                        let rc = node_id[..node_id.len()-1].to_string() + "-";
                        graph.nodes.remove(rc.as_str());
                    } else if node_id.ends_with('-') {
                        let rc = node_id[..node_id.len()-1].to_string() + "+";
                        graph.nodes.remove(rc.as_str());
                    }
                }
                println!("Removed {} oriented nodes from small components (<2)", small_comp_count);

                // tip trimming
                let before_trim = graph.nodes.len();
                tip_trimming::trim_tips(&mut graph, max_tip_len);
                let after_trim = graph.nodes.len();
                println!("Removed {} nodes by tip trimming", before_trim.saturating_sub(after_trim));

            }

            // heuristic simplification
            println!("Applying heuristic simplification: removing weak edges...");
            heuristic_simplification::remove_weak(&mut graph);

            // write graph snapshot into output dir
            let dot_path = out_dir.join("overlap_before_compression.dot");
            let dot_str = dot_path.to_str().ok_or("invalid output path")?;
            graph.write_dot(dot_str)?;

            // compress into unitigs into output dir
            let out_path = out_dir.join(format!("{}.fa", args.output_prefix));
            let out_str = out_path.to_str().ok_or("invalid output path")?;
            let compressed = compress_graph::compress_unitigs(&graph, &args.reads_fq, out_str);
            println!("Assembly produced {} unitigs (written to {})", compressed.unitigs.len(), out_str);
            let gfa_path = out_dir.join(format!("{}.gfa", args.output_prefix));
            let gfa_str = gfa_path.to_str().ok_or("invalid output path")?;
            compressed.write_gfa(gfa_str)?;
            println!("Wrote GFA to {}", gfa_str);
        }
    }

    Ok(())
}