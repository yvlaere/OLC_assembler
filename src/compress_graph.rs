/// graph compression module
/// creates a compressed graph of unitigs from an overlap graph
/// 1. get the indegree and outdegree of each node
/// 2. get non-circular unitigs (start at nodes with indegree != 1 or outdegree != 1)
/// 3. get circular unitigs (remaining unvisited nodes)

use std::collections::{HashMap, HashSet};

pub struct UnitigMember {
    pub node_id: String,
    // node id and overlap length to the next node in the unitig
    pub edge: (String, u32),
}

pub struct Unitig {
    pub id: usize,
    pub members: Vec<UnitigMember>,
    pub cached_seq: Option<String>,
}

pub struct UnitigEdge {
    pub from: usize,
    pub to: usize,
    pub overlap: u32,
}

pub struct CompressedGraph {
    pub unitigs: Vec<Unitig>,
    //pub edges: Vec<UnitigEdge>,
}

/// Main function: compress maximal non-branching paths into unitigs.
/// Preserves member lists and the overlap lengths between them.
pub fn compress_unitigs(graph: &crate::create_overlap_graph::OverlapGraph,) -> CompressedGraph {
    
    // 1) create a map of indegrees
    let mut indegree: HashMap<String, usize> = HashMap::new();
    for id in graph.nodes.keys() {
        indegree.insert(id.clone(), 0);
    }
    for (source_id, node) in &graph.nodes {
        for e in &node.edges {
            *indegree.entry(e.target_id.clone()).or_default() += 1;
        }
    }

    let mut visited: HashSet<String> = HashSet::new();
    let mut unitigs: Vec<Unitig> = Vec::new();

    // Helper to extract the single outgoing neighbor if outdeg == 1
    let out_single = |g: &crate::create_overlap_graph::OverlapGraph, cur: &str| -> Option<(String, u32)> {
        g.nodes.get(cur).and_then(|n| {
            if n.edges.len() == 1 {
                let e = &n.edges[0];
                Some((e.target_id.clone(), e.edge_len))
            } 
            else {
                None
            }
        })
    };

    // 2) non-circular unitigs, start unitigs at nodes where indegree != 1 || outdeg != 1
    for (id, node) in &graph.nodes {
        let indegree_i = *indegree.get(id).unwrap_or(&0);
        let outdeg_i = node.edges.len();

        // skip if already visited
        if visited.contains(id) { println!("already visited"); continue; }

        // start a unitig if indegree != 1 or outdegree != 1 (i.e., not a simple internal node)
        if indegree_i != 1 || outdeg_i != 1 {

            println!("starting unitig at node {}", id);
            println!("indegree: {}, outdegree: {}", indegree_i, outdeg_i);

            // create a unitig for each outgoing edge
            // if outdeg == 0, there will be no unitig
            for out_edge_i in 0..outdeg_i {

                println!("processing outgoing edge {}", out_edge_i);

                // start a new unitig from id
                let mut members: Vec<UnitigMember> = Vec::new();
                let mut cur = id.clone();
                //members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
                visited.insert(cur.clone());

                // check the next outgoing edge
                let (second, overlap) =  {
                    let e = &node.edges[out_edge_i];
                    (e.target_id.clone(), e.edge_len)
                };

                println!("checking second node {}", second);
                let second_indegree = *indegree.get(&second).unwrap_or(&0);
                // don't add the node that breaks the chain
                if second_indegree != 1 { println!("second node indegree != 1"); continue; }
                // stop if second is already visited
                if visited.contains(&second) { println!("second node already visited"); continue; }
                // push the first node into the unitig members
                members.push(UnitigMember { node_id: cur.clone(), edge: (second.clone(), overlap) });
                visited.insert(second.clone());
                cur = second.clone();
                
                // extend forward from second untill the end
                while let Some((next, overlap)) = out_single(graph, &cur) {
                    println!("extending unitig");
                    let next_indegree = *indegree.get(&next).unwrap_or(&0);
                    // don't add the node that breaks the chain
                    if next_indegree != 1 { break; }
                    // stop if next is already visited
                    if visited.contains(&next) { break; }
                    // push cur to the unitig members
                    members.push(UnitigMember { node_id: cur.clone(), edge: (next.clone(), overlap) });
                    visited.insert(next.clone());
                    cur = next;
                }

                // create the unitig
                let uid = unitigs.len();
                unitigs.push(Unitig { id: uid, members, cached_seq: None });
            }
        }
    }

    println!("\n");

    // 3) circular unitigs, handle remaining nodes that are still unvisited
    for id in graph.nodes.keys() {
        if visited.contains(id) { continue; }
        println!("starting circular unitig at node {}", id);

        // start a circular unitig
        let mut cur = id.clone();
        let mut members: Vec<UnitigMember> = Vec::new();

        loop {
            visited.insert(cur.clone());
            // to follow, get the unique outgoing edge
            println!("extending circular unitig");
            let next_edge = out_single(graph, &cur);
            if next_edge.is_none() {
                // shouldn't happen in pure cycle, break defensively
                // members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
                println!("warning: broken cycle detected during unitig compression");
                break;
            }
            let (next, overlap) = next_edge.unwrap();

            let next_indegree = *indegree.get(&next).unwrap_or(&0);
            // don't add the node that breaks the chain
            if next_indegree != 1 { break; }
            // stop if next is already visited
            if visited.contains(&next) { break; }
            // push cur to the unitig members
            members.push(UnitigMember { node_id: cur.clone(), edge: (next.clone(), overlap) });

            cur = next;
        }

        // register unitig
        let uid = unitigs.len();
        unitigs.push(Unitig { id: uid, members, cached_seq: None });
    }

    CompressedGraph { unitigs }
}