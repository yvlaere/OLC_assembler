/// graph compression module
/// creates a compressed graph of unitigs from an overlap graph
/// 1. get the indegree and outdegree of each node
/// 2. get non-circular unitigs (start at nodes with indegree != 1 or outdegree != 1)
/// 3. get circular unitigs (remaining unvisited nodes)

use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use crate::utils;

pub struct UnitigMember {
    pub node_id: String,
    // id of target node and edge length to that node
    pub edge: (String, u32),
}

pub struct Unitig {
    pub id: usize,
    pub members: Vec<UnitigMember>,
    pub fasta_seq: Option<String>,
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
pub fn compress_unitigs(graph: &crate::create_overlap_graph::OverlapGraph, fastq_path: &str, fasta_path: &str) -> CompressedGraph {
    
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
                let (second, edge_len) =  {
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
                members.push(UnitigMember { node_id: cur.clone(), edge: (second.clone(), edge_len) });
                visited.insert(second.clone());
                cur = second.clone();
                
                // extend forward from second untill the end
                while let Some((next, edge_len)) = out_single(graph, &cur) {
                    println!("extending unitig");
                    let next_indegree = *indegree.get(&next).unwrap_or(&0);
                    // don't add the node that breaks the chain
                    if next_indegree != 1 { break; }
                    // stop if next is already visited
                    if visited.contains(&next) { break; }
                    // push cur to the unitig members
                    members.push(UnitigMember { node_id: cur.clone(), edge: (next.clone(), edge_len) });
                    visited.insert(next.clone());
                    cur = next;
                }

                // create the unitig
                let uid = unitigs.len();
                unitigs.push(Unitig { id: uid, members, fasta_seq: None });
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
                println!("warning: broken cycle detected during unitig compression");
                break;
            }
            let (next, edge_len) = next_edge.unwrap();

            let next_indegree = *indegree.get(&next).unwrap_or(&0);
            // don't add the node that breaks the chain
            if next_indegree != 1 { break; }
            // stop if next is already visited
            if visited.contains(&next) { break; }
            // push cur to the unitig members
            members.push(UnitigMember { node_id: cur.clone(), edge: (next.clone(), edge_len) });

            cur = next;
        }

        // register unitig
        let uid = unitigs.len();
        unitigs.push(Unitig { id: uid, members, fasta_seq: None });
    }

    // load fastq sequences
    println!("Loading FASTQ sequences from {}...", fastq_path);
    let fastq_seqs = load_fastq_sequences(fastq_path).unwrap();

    // generate fasta sequences for unitigs and write to fasta_path
    println!("Generating unitig sequences and writing to FASTA...");
    for unitig in unitigs.iter_mut() {
        let seq = unitig_sequence(unitig, graph, fastq_seqs.clone()).unwrap();
        unitig.fasta_seq = Some(seq);
    }
    // write to fasta file
    {
        let mut fasta_file = std::fs::File::create(fasta_path).unwrap();
        for unitig in unitigs.iter() {
            let header = format!(">unitig_{} len={}\n", unitig.id, unitig.members.len());
            let seq = unitig.fasta_seq.as_ref().unwrap();
            use std::io::Write;
            fasta_file.write_all(header.as_bytes()).unwrap();
            fasta_file.write_all(seq.as_bytes()).unwrap();
            fasta_file.write_all(b"\n").unwrap();
        }
    }

    CompressedGraph { unitigs }
}

fn load_fastq_sequences(fastq_path: &str) -> Result<HashMap<String, String>, String> {
    let mut seq_map: HashMap<String, String> = HashMap::new();

    let reader = match std::fs::File::open(fastq_path) {
        Ok(f) => f,
        Err(e) => return Err(format!("failed to open FASTQ file '{}': {}", fastq_path, e)),
    };
    let buf_reader = std::io::BufReader::new(reader);
    let mut lines = buf_reader.lines();

    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('@') {
            return Err(format!("invalid FASTQ format: expected header line starting with '@', got '{}'", header));
        }
        let seq = match lines.next() {
            Some(Ok(s)) => s,
            _ => return Err("invalid FASTQ format: missing sequence line".to_string()),
        };
        // skip plus line
        lines.next();
        // skip quality line
        lines.next();

        let id = header[1..].split_whitespace().next().unwrap().to_string();
        seq_map.insert(id, seq);
    }

    Ok(seq_map)
}

/// Build nucleotide sequence for `unitig` by concatenating node sequences and
/// removing overlaps recorded in UnitigMember.edge.(target, edge_len).
pub fn unitig_sequence(unitig: &Unitig, graph: &crate::create_overlap_graph::OverlapGraph, fastq_seqs: HashMap<String, String>) -> Result<String, String> {
    if unitig.members.is_empty() {
       println!("unitig has no members; cannot infer sequence");
    }

    // Helper to get sequence for a node id
        let get_seq = |node_id: &str| -> Result<String, String> {
            // Node ID must end with '+' or '-'
            let orientation = node_id.chars().last().ok_or_else(|| {
                format!("node id '{}' has no orientation suffix (+/-)", node_id)
            })?;
            let read_id = &node_id[..node_id.len() - 1];
            let node = graph.nodes.get(node_id)
                .ok_or_else(|| format!("node_id '{}' not found in overlap graph", node_id))?;
            let seq = fastq_seqs.get(read_id)
                .ok_or_else(|| format!("sequence for read_id '{}' not found in FASTQ sequences", read_id))?;
            match orientation {
                '+' => Ok(seq.clone()),
                '-' => Ok(utils::rev_comp(seq)),
                _ => Err(format!("invalid orientation '{}' in node id '{}'", orientation, node_id)),
            }
        };

    let mut out = String::new();

    // For each member, append the prefix untill the target (skip overlap bases)
    for member in &unitig.members {
        // an edge describes the target node and the edge length of the original node to reach that target
        let (target_id, edge_len) = &member.edge;
        let seq = get_seq(&member.node_id)?;
        let edge_len_usize = *edge_len as usize;

        if edge_len_usize > seq.len() {
            return Err(format!("edge length ({}) larger than seq length ({}) for node {} to target {}", edge_len_usize, seq.len(), member.node_id, target_id));
        }

        // Append the non-overlapping prefix
        out.push_str(&seq[..edge_len_usize]);
    }

    Ok(out)
}