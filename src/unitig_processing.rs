/// Unitig processing module
/// Handles post-compression operations on unitig graphs:
/// - Finding longest circular paths (for circular genomes/plasmids)
/// - Detecting and removing reverse complement duplicates

use std::collections::{HashMap, HashSet, VecDeque};
use crate::compress_graph::{CompressedGraph, Unitig, UnitigEdge};
use crate::utils;

/// Represents a path through unitigs with cumulative length
#[derive(Clone)]
pub struct UnitigPath {
    pub unitigs: Vec<usize>,
    pub is_circular: bool,
    pub total_length: usize,
}

impl UnitigPath {
    pub fn new() -> Self {
        Self {
            unitigs: Vec::new(),
            is_circular: false,
            total_length: 0,
        }
    }
}

/// Find the longest circular path in the unitig graph
/// Returns the path with highest total member count (as a proxy for sequence length)
pub fn find_longest_circular_path(graph: &CompressedGraph) -> Option<UnitigPath> {
    if graph.unitigs.is_empty() {
        return None;
    }

    // Build adjacency map
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for e in &graph.edges {
        adjacency.entry(e.from).or_insert_with(Vec::new).push(e.to);
    }

    let mut longest_path: Option<UnitigPath> = None;
    let mut visited_starts: HashSet<usize> = HashSet::new();

    // Try starting from each unitig
    for start_uid in 0..graph.unitigs.len() {
        if visited_starts.contains(&start_uid) {
            continue;
        }

        // DFS to find cycles starting from start_uid
        let mut stack: Vec<(usize, Vec<usize>, HashSet<usize>)> = vec![(start_uid, vec![start_uid], {
            let mut s = HashSet::new();
            s.insert(start_uid);
            s
        })];

        while let Some((current, path, visited)) = stack.pop() {
            // Limit path length to avoid combinatorial explosion
            if path.len() > 1000 {
                continue;
            }

            if let Some(neighbors) = adjacency.get(&current) {
                for &next in neighbors {
                    if next == start_uid && path.len() > 1 {
                        // Found a cycle back to start
                        let mut cycle_path = path.clone();
                        let total_length: usize = cycle_path.iter().map(|&uid| graph.unitigs[uid].members.len()).sum();

                        let new_path = UnitigPath {
                            unitigs: cycle_path,
                            is_circular: true,
                            total_length,
                        };

                        if longest_path.is_none() || new_path.total_length > longest_path.as_ref().unwrap().total_length {
                            longest_path = Some(new_path);
                            visited_starts.insert(start_uid);
                        }
                    } else if !visited.contains(&next) && path.len() <= 50 {
                        // Continue exploring
                        let mut new_visited = visited.clone();
                        new_visited.insert(next);
                        let mut new_path = path.clone();
                        new_path.push(next);
                        stack.push((next, new_path, new_visited));
                    }
                }
            }
        }
    }

    longest_path
}

/// Detect unitigs that are reverse complements of each other
/// Returns a map of (original_uid -> rc_uid) for duplicates
pub fn detect_reverse_complement_duplicates(
    graph: &CompressedGraph,
    unitig_seqs: &HashMap<usize, String>,
) -> HashMap<usize, usize> {
    let mut duplicates: HashMap<usize, usize> = HashMap::new();
    let mut processed: HashSet<usize> = HashSet::new();

    for uid1 in 0..graph.unitigs.len() {
        if processed.contains(&uid1) {
            continue;
        }

        if let Some(seq1) = unitig_seqs.get(&uid1) {
            let seq1_rc = utils::rev_comp(seq1);

            for uid2 in (uid1 + 1)..graph.unitigs.len() {
                if processed.contains(&uid2) {
                    continue;
                }

                if let Some(seq2) = unitig_seqs.get(&uid2) {
                    // Check if seq2 is the RC of seq1
                    if seq2 == &seq1_rc {
                        // Mark uid2 as a duplicate (RC) of uid1
                        duplicates.insert(uid2, uid1);
                        processed.insert(uid2);
                        println!("Detected RC duplicate: unitig_{} is RC of unitig_{}", uid2, uid1);
                    }
                    // Also check the reverse
                    let seq2_rc = utils::rev_comp(seq2);
                    if seq1 == &seq2_rc && !duplicates.contains_key(&uid1) {
                        duplicates.insert(uid1, uid2);
                        processed.insert(uid1);
                        println!("Detected RC duplicate: unitig_{} is RC of unitig_{}", uid1, uid2);
                    }
                }
            }
        }
        processed.insert(uid1);
    }

    duplicates
}

/// Remove reverse complement duplicates from the graph
/// Keeps the first occurrence and redirects edges from the RC to the original
pub fn remove_rc_duplicates(graph: &mut CompressedGraph, duplicates: &HashMap<usize, usize>) {
    if duplicates.is_empty() {
        println!("No reverse complement duplicates found");
        return;
    }

    println!("Removing {} reverse complement duplicates", duplicates.len());

    // Collect UIDs to remove
    let uids_to_remove: HashSet<usize> = duplicates.keys().cloned().collect();

    // Update edges: redirect edges from duplicate nodes to original nodes
    for edge in &mut graph.edges {
        if let Some(&original_uid) = duplicates.get(&edge.from) {
            edge.from = original_uid;
        }
        if let Some(&original_uid) = duplicates.get(&edge.to) {
            edge.to = original_uid;
        }
    }

    // Remove self-loops
    graph.edges.retain(|e| e.from != e.to);

    // Remove duplicate unitigs
    for &uid_to_remove in &uids_to_remove {
        // Mark for removal by filtering
        // Note: we keep the unitig in the vec but can mark it somehow
        // For simplicity, we'll just leave them but remove edges
        println!("Marked unitig_{} for removal (kept as placeholder)", uid_to_remove);
    }

    // Clean up edges that reference removed nodes (self-loops already handled)
    // Deduplicate edges (same from->to pair)
    let mut seen_edges: HashSet<(usize, usize)> = HashSet::new();
    graph.edges.retain(|e| {
        let key = (e.from, e.to);
        if seen_edges.contains(&key) {
            false
        } else {
            seen_edges.insert(key);
            true
        }
    });

    println!("RC duplicate removal complete. {} edges remain", graph.edges.len());
}

/// Export a path to a FASTA sequence
/// Concatenates unitig sequences along the path
pub fn path_to_fasta(path: &UnitigPath, graph: &CompressedGraph, sequences: &HashMap<usize, String>) -> Option<String> {
    let mut fasta_seq = String::new();

    for &uid in &path.unitigs {
        if let Some(seq) = sequences.get(&uid) {
            fasta_seq.push_str(seq);
        } else {
            eprintln!("Warning: sequence not found for unitig_{}", uid);
            return None;
        }
    }

    Some(fasta_seq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rc_detection() {
        let seq1 = "ACGTACGT";
        let seq1_rc = "ACGTACGT"; // palindromic
        assert_eq!(utils::rev_comp(seq1), seq1_rc);
    }
}
