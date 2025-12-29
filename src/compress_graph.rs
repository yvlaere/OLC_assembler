/// graph compression module
/// creates a compressed graph of unitigs from an overlap graph

use std::collections::{HashMap, HashSet};

pub struct UnitigMember {
    pub node_id: String,
    pub overlap_from_prev: u32,
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
    pub edges: Vec<UnitigEdge>,
    pub node_to_unitig: HashMap<String, usize>,
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
    let mut node_to_unitig: HashMap<String, usize> = HashMap::new();

    // Helper to extract the single outgoing neighbor if outdeg == 1
    let out_single = |g: &crate::create_overlap_graph::OverlapGraph, cur: &str| -> Option<(String, u32)> {
        g.nodes.get(cur).and_then(|n| {
            if n.edges.len() == 1 {
                let e = &n.edges[0];
                Some((e.target_id.clone(), e.edge_len))
            } else {
                None
            }
        })
    };

    // 2) Start unitigs at nodes where indegree != 1 || outdeg != 1
    for (id, node) in &graph.nodes {
        let indegree_i = *indegree.get(id).unwrap_or(&0);
        let outdeg_i = node.edges.len();

        if indegree_i != 1 || outdeg_i != 1 {
            if visited.contains(id) { continue; }
            // start a new unitig from id
            println!("starting unitig at node {}", id);
            let mut members: Vec<UnitigMember> = Vec::new();
            let mut cur = id.clone();
            members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
            visited.insert(cur.clone());

            // walk forward while next node has indegree == 1 && outdeg == 1
            while let Some((next, overlap)) = out_single(graph, &cur) {
                println!("extending unitig");
                let next_indegree = *indegree.get(&next).unwrap_or(&0);
                if next_indegree != 1 {
                    // still add the next as the last node? no — stop here.
                    break;
                }
                if visited.contains(&next) { break; }
                // push next with overlap_from_prev = overlap (cur -> next)
                members.push(UnitigMember { node_id: next.clone(), overlap_from_prev: overlap });
                visited.insert(next.clone());
                cur = next;
            }

            let uid = unitigs.len();
            for m in &members {
                node_to_unitig.insert(m.node_id.clone(), uid);
            }
            unitigs.push(Unitig { id: uid, members, cached_seq: None });
        }
    }

    // 3) Handle remaining nodes that are still unvisited (they belong to cycles where indegree==1 and outdeg==1)
    for id in graph.nodes.keys() {
        if visited.contains(id) { continue; }
        // start a cycle unitig
        let mut cur = id.clone();
        let mut members: Vec<UnitigMember> = Vec::new();

        loop {
            visited.insert(cur.clone());
            // to follow, get the unique outgoing edge
            let next_edge = out_single(graph, &cur);
            if next_edge.is_none() {
                // shouldn't happen in pure cycle, break defensively
                members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
                break;
            }
            let (next, overlap) = next_edge.unwrap();

            // add current with proper overlap_from_prev:
            if members.is_empty() {
                members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
            } else {
                // calculate overlap_from_prev for cur: we need the overlap from previous->cur
                // but in the cycle walk we added that when we pushed the node previously
            }

            // move to next; stop when next is already in members (closed the cycle)
            if visited.contains(&next) {
                // we've reached a previously visited node — close cycle
                // push next only if not included
                if !members.iter().any(|m| m.node_id == next) {
                    // find overlap from cur->next
                    members.push(UnitigMember { node_id: next.clone(), overlap_from_prev: overlap });
                }
                break;
            }

            // push next with overlap_from_prev = overlap (cur->next)
            members.push(UnitigMember { node_id: next.clone(), overlap_from_prev: overlap });
            cur = next;
        }

        // register unitig (we inserted members maybe with duplicated last node in cycle; that's acceptable)
        let uid = unitigs.len();
        for m in &members {
            node_to_unitig.insert(m.node_id.clone(), uid);
        }
        unitigs.push(Unitig { id: uid, members, cached_seq: None });
    }

    // 4) Build unitig-to-unitig edges by mapping outgoing edges of last member of each unitig
    let mut edges: Vec<UnitigEdge> = Vec::new();
    let mut seen_pairs: HashSet<(usize, usize, u32)> = HashSet::new(); // dedupe exact triples

    for u in &unitigs {
        if u.members.is_empty() { continue; }
        let last = u.members.last().unwrap();
        if let Some(last_node) = graph.nodes.get(&last.node_id) {
            for e in &last_node.edges {
                let tgt_node = &e.target_id;
                let overlap = e.edge_len;
                if let Some(&tgt_uid) = node_to_unitig.get(tgt_node) {
                    if tgt_uid == u.id { continue; } // skip self-loops if you want to ignore them
                    let triple = (u.id, tgt_uid, overlap);
                    if !seen_pairs.contains(&triple) {
                        edges.push(UnitigEdge { from: u.id, to: tgt_uid, overlap: overlap });
                        seen_pairs.insert(triple);
                    }
                }
            }
        }
    }

    CompressedGraph { unitigs, edges, node_to_unitig }
}

/// Build the sequence for a unitig on demand using the original read sequences and overlaps.
/// - read_seqs: map "node_id" -> seq (if node_id encodes orientation like "read1+"/ "read1-",
///   the caller should provide oriented sequences or we detect '-' suffix and RC).
pub fn build_unitig_sequence(
    unitig: &mut Unitig,
    graph: &crate::create_overlap_graph::OverlapGraph,
    read_seqs: &HashMap<String, String>,
    force_recompute: bool,
) -> Option<String> {
    if !force_recompute {
        if let Some(seq) = &unitig.cached_seq {
            return Some(seq.clone());
        }
    }

    if unitig.members.is_empty() { unitig.cached_seq = Some(String::new()); return Some(String::new()); }

    // helper RC if you encode orientation in node id (e.g., "name-")
    fn rc_if_needed(node_id: &str, seq: &str) -> String {
        if node_id.ends_with('-') {
            reverse_complement(seq)
        } else {
            seq.to_owned()
        }
    }

    let first = &unitig.members[0];
    let first_seq = read_seqs.get(&first.node_id)?;
    let mut assembled = rc_if_needed(&first.node_id, first_seq);

    for m in unitig.members.iter().skip(1) {
        let cur_seq_raw = read_seqs.get(&m.node_id)?;
        let cur_seq = rc_if_needed(&m.node_id, cur_seq_raw);
        let ov = m.overlap_from_prev as usize;
        if ov > cur_seq.len() {
            // inconsistent overlap: abort
            return None;
        }
        // append non-overlapping suffix
        assembled.push_str(&cur_seq[ov..]);
    }

    unitig.cached_seq = Some(assembled.clone());
    Some(assembled)
}

/// Reverse complement for ASCII DNA sequences (A,C,G,T,N)
fn reverse_complement(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.as_bytes().iter().rev() {
        let rc = match *c as char {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            'N' => 'N',
            'n' => 'n',
            other => other, // passthrough
        };
        out.push(rc);
    }
    out
}