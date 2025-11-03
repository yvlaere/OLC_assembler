use std::collections::{HashMap, HashSet};

/// Velvet-style tip trimming for an OverlapGraph.
/// - max_tip_nodes: maximum number of nodes a tip may contain to be eligible (use 4).
/// Returns the number of nodes removed.

// recompute indegree, outdegree, edge multiplicity, and reverse adjacency.
fn recompute_stats(graph: &crate::create_overlap_graph::OverlapGraph,) -> (
    HashMap<String, usize>,                         // indeg
    HashMap<String, usize>,                         // outdeg
    HashMap<(String, String), usize>,               // edge multiplicity (src, tgt) -> count
    HashMap<String, Vec<String>>,                   // predecessors: tgt -> Vec<src>
) {
    let mut indeg: HashMap<String, usize> = HashMap::new();
    let mut outdeg: HashMap<String, usize> = HashMap::new();
    let mut edge_mult: HashMap<(String, String), usize> = HashMap::new();
    let mut preds: HashMap<String, Vec<String>> = HashMap::new();

    // initialize indeg/outdeg to zero for all nodes
    for id in graph.nodes.keys() {
        indeg.insert(id.clone(), 0);
        outdeg.insert(id.clone(), 0);
        preds.insert(id.clone(), Vec::new());
    }

    for (src, node) in graph.nodes.iter() {
        let o = node.edges.len();
        outdeg.insert(src.clone(), o);
        for (tgt, _len) in node.edges.iter() {
            // increment indeg
            *indeg.entry(tgt.clone()).or_default() += 1;
            // edge multiplicity: count each (src,tgt) occurrence
            *edge_mult.entry((src.clone(), tgt.clone())).or_default() += 1;
            // reverse adjacency
            preds.entry(tgt.clone()).or_default().push(src.clone());
        }
    }

    (indeg, outdeg, edge_mult, preds)
}

// remove a set of nodes (by id) from the graph and remove references to them in other nodes.
// compute reverse-complement node id for oriented node ids ending in '+' or '-'
fn rc_node(id: &str) -> String {
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

// remove a set of nodes (by id) from the graph and remove references to them in other nodes.
// Also remove reverse-complement counterparts to keep the bigraph synchronized.
// Returns the number of nodes actually removed (including RC counterparts).
fn remove_nodes(graph: &mut crate::create_overlap_graph::OverlapGraph, to_remove: &HashSet<String>) -> usize {
    // expand set to include reverse complements
    let mut expanded: HashSet<String> = HashSet::new();
    for id in to_remove.iter() {
        expanded.insert(id.clone());
        expanded.insert(rc_node(id));
    }

    // Remove nodes from map
    for id in expanded.iter() {
        graph.nodes.remove(id);
    }
    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges.retain(|(tgt, _len)| !expanded.contains(tgt));
    }

    expanded.len()
}

pub fn trim_tips(graph: &mut crate::create_overlap_graph::OverlapGraph, max_tip_nodes: usize,) -> usize {
    
    // No tip trimming
    if max_tip_nodes == 0 {
        return 0;
    }

    println!("Starting tip trimming with max_tip_nodes={}", max_tip_nodes);

    let mut counter = 0;

    // Main iterative loop
    let mut total_removed = 0usize;

    loop {
        let (indeg, outdeg, edge_mult, preds) = recompute_stats(graph);

        println!("Counter {}: Graph has {} nodes", counter, graph.nodes.len());
        counter += 1;

        // Build a list of current tip nodes: nodes with outdeg == 0 OR indeg == 0
        // We'll create candidates for both "tail tips" (outdeg==0) and "head tips" (indeg==0).
        let mut tip_candidates: Vec<(String, bool)> = Vec::new(); // (node_id, is_tail_tip)
        for id in graph.nodes.keys() {
            let id_s = id.clone();
            if *outdeg.get(id).unwrap_or(&0) == 0 {
                tip_candidates.push((id_s.clone(), true)); // tail tip (dead end at out)
            } else if *indeg.get(id).unwrap_or(&0) == 0 {
                tip_candidates.push((id_s.clone(), false)); // head tip (dead end at in)
            }
        }

        if tip_candidates.is_empty() {
            println!("No tip candidates found, stopping.");
            break;
        }

        // For each candidate tip, attempt to trace the maximal non-branching chain up to max_tip_nodes
        // and record the multiplicity of the junction arc (arc from junction -> first tip node for tail tips,
        // or arc from first tip node -> junction for head tips).
        #[derive(Clone)]
        struct TipInfo {
            tip_node: String,
            is_tail: bool,
            chain: Vec<String>,     // ordered from junction-neighbor ... tip_node, or from tip -> junction neighbor depending
            junction: String,
            junction_arc_mult: usize,
        }

        let mut candidates_info: Vec<TipInfo> = Vec::new();

        for (tip_node, is_tail) in tip_candidates.into_iter() {
            // Walk outward from tip_node toward junction.
            // For tail tips (outdeg==0) we walk backwards following unique predecessor while those nodes are non-branching.
            // For head tips (indeg==0) we walk forwards following unique successor while non-branching.
            let mut chain: Vec<String> = Vec::new();
            chain.push(tip_node.clone());

            if is_tail {
                // walk backward: while the current node has exactly one predecessor p and that p has outdeg==1 and indeg==1
                let mut cur = tip_node.clone();
                for _ in 0..(max_tip_nodes - 1) { // -1 because we already have the tip node in chain
                    // get predecessors of cur
                    let pred_list = preds.get(&cur).map(|v| v.clone()).unwrap_or_else(Vec::new);
                    if pred_list.len() != 1 {
                        break;
                    }
                    let p = pred_list[0].clone();
                    let p_out = *outdeg.get(&p).unwrap_or(&0);
                    let p_in = *indeg.get(&p).unwrap_or(&0);
                    // continue only if p is non-branching (in==1 and out==1)
                    if p_out == 1 && p_in == 1 {
                        chain.push(p.clone());
                        cur = p;
                    } else {
                        // p is the junction (it may have other outgoing arcs)
                        // include p in chain? Velvet considers chain up to node adjacent to junction; here we stop before adding p,
                        // but we need to identify junction node separately.
                        break;
                    }
                }
                // determine junction node: node adjacent to chain that connects to rest of graph
                // For tail tip, the junction is the unique predecessor of the last node in chain (if exists),
                // or the last node in chain's predecessor that is branching. We treat as:
                let last = chain.last().unwrap().clone();
                // find predecessor(s) of last
                let pred_list = preds.get(&last).map(|v| v.clone()).unwrap_or_else(Vec::new);
                // choose a predecessor that is outside the chain (the junction). If multiple exist, pick the first.
                let junction_opt = pred_list.into_iter().find(|p| !chain.contains(p));
                let junction = match junction_opt {
                    Some(p) => p,
                    None => continue, // no external predecessor found -> skip
                };
                // compute multiplicity of the arc junction -> neighbor_on_chain
                // neighbor_on_chain is the node in chain that is adjacent to junction; that is the previous element if chain contains it
                // we need to find which element in chain is connected from junction
                let neighbor_on_chain = {
                    // The neighbor should be the node in chain that has junction as predecessor; usually the last element of chain.
                    chain.last().unwrap().clone()
                };
                let arc_mult = *edge_mult.get(&(junction.clone(), neighbor_on_chain.clone())).unwrap_or(&0usize);

                // also ensure chain length is <= max_tip_nodes
                if chain.len() <= max_tip_nodes {
                    candidates_info.push(TipInfo {
                        tip_node: tip_node.clone(),
                        is_tail: true,
                        chain: chain.clone(),
                        junction: junction.clone(),
                        junction_arc_mult: arc_mult,
                    });
                }
            } else {
                // head tip: indeg == 0. Walk forward following unique outgoing edges while nodes are non-branching.
                let mut cur = tip_node.clone();
                for _ in 0..(max_tip_nodes - 1) {
                    // get successors of cur
                    let succs: Vec<String> = graph.nodes.get(&cur)
                        .map(|n| n.edges.iter().map(|(t, _)| t.clone()).collect())
                        .unwrap_or_else(Vec::new);
                    if succs.len() != 1 {
                        break;
                    }
                    let s = succs[0].clone();
                    let s_out = *outdeg.get(&s).unwrap_or(&0);
                    let s_in = *indeg.get(&s).unwrap_or(&0);
                    if s_out == 1 && s_in == 1 {
                        chain.push(s.clone());
                        cur = s;
                    } else {
                        // s is the junction (branching)
                        break;
                    }
                }
                // determine junction node: the unique successor of the last chain node
                let last = chain.last().unwrap().clone();
                let succs = graph.nodes.get(&last)
                    .map(|n| n.edges.iter().map(|(t,_len)| t.clone()).collect::<Vec<_>>())
                    .unwrap_or_else(Vec::new);
                // pick a successor that is outside the chain (junction). If none, skip.
                let junction_opt = succs.into_iter().find(|s| !chain.contains(s));
                let junction = match junction_opt {
                    Some(s) => s,
                    None => continue,
                };
                // arc multiplicity: arc from neighbor_on_chain -> junction (neighbor_on_chain is last element of chain)
                let neighbor_on_chain = last.clone();
                let arc_mult = *edge_mult.get(&(neighbor_on_chain.clone(), junction.clone())).unwrap_or(&0usize);

                if chain.len() <= max_tip_nodes {
                    candidates_info.push(TipInfo {
                        tip_node: tip_node.clone(),
                        is_tail: false,
                        chain: chain.clone(),
                        junction: junction.clone(),
                        junction_arc_mult: arc_mult,
                    });
                }
            }
        } // end building candidates_info

        if candidates_info.is_empty() {
            println!("No valid tip candidates after tracing, stopping.");
            break;
        }

        println!("Found {} tip candidates for removal consideration", candidates_info.len());

        // Previously we applied a minority-count / length check at junctions to select eligible tips.
        // That test caused many legitimate tips to be kept; remove that complexity and treat all traced
        // candidates as eligible (subject to basic sanity checks later).
        let mut eligible: Vec<TipInfo> = candidates_info.into_iter()
            .filter(|info| {
                // require the junction node still exists and has at least one outgoing arc
                graph.nodes.get(&info.junction).map(|jn| !jn.edges.is_empty()).unwrap_or(false)
            })
            .collect();

        if eligible.is_empty() {
            println!("No tip candidates eligible for removal, stopping.");
            break;
        }

        println!("{} tip candidates eligible for removal", eligible.len());

        // Sort eligible by junction_arc_mult ascending (Velvet removes low multiplicity first)
        eligible.sort_by_key(|x| x.junction_arc_mult);

        // We'll attempt to remove candidates in order, removing the chain nodes for each eligible candidate
        // But we must ensure we don't remove nodes that are no longer present or have become non-tips due to earlier removals in this iteration.
        let mut removed_in_this_round: HashSet<String> = HashSet::new();
        for info in eligible.into_iter() {
            // sanity checks: chain nodes must still exist in graph and still be tips (outdeg==0 or indeg==0)
            let mut all_present = true;
            for n in info.chain.iter() {
                if !graph.nodes.contains_key(n) {
                    all_present = false;
                    break;
                }
            }
            if !all_present {
                continue;
            }

            // re-check that the chain still forms a proper tip:
            // tail-tip: last node in chain should have outdeg==0
            // head-tip: first node in chain should have indeg==0
            if info.is_tail {
                let tip_last = info.tip_node.clone();
                // tip node must still have outdeg == 0
                let outd = graph.nodes.get(&tip_last).map(|n| n.edges.len()).unwrap_or(0);
                if outd != 0 {
                    continue;
                }
            } else {
                let tip_first = info.tip_node.clone();
                // tip node must still have indeg == 0
                let indeg_map = {
                    // compute indegree for the tip_first quickly by scanning predecessors (cheap for single check)
                    let mut c = 0usize;
                    for (_src, n) in graph.nodes.iter() {
                        for (tgt, _len) in n.edges.iter() {
                            if tgt == &tip_first {
                                c += 1;
                            }
                        }
                    }
                    c
                };
                if indeg_map != 0 {
                    continue;
                }
            }

            // Remove all chain nodes (chain contains the tip node and possibly up-to max_tip_nodes-1 upstream/downstream nodes)
            for n in info.chain.iter() {
                removed_in_this_round.insert(n.clone());
            }

            // Update total_removed lazily after batch removal
        }

        if removed_in_this_round.is_empty() {
            // no actual removals possible
            break;
        }

    // Perform actual removal: delete nodes and purge incoming edges
    let batch = removed_in_this_round.clone();
    let removed_count = remove_nodes(graph, &batch);
    total_removed += removed_count;

        // After removals, loop again (recompute stats in next iteration)
    } // end main loop

    total_removed
}