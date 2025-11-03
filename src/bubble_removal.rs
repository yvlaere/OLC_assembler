use std::collections::{HashMap, HashSet, VecDeque};

/// Bubble removal using a "tour bus" style search.
///
/// This module provides a conservative bubble removal routine suitable for
/// overlap graphs represented by `crate::create_overlap_graph::OverlapGraph`.
///
/// Public API:
/// - remove_bubbles(graph, max_bubble_len, min_support_ratio) -> usize
///
/// Behavior summary:
/// - For each node `u` with out-degree >= 2, consider every pair of distinct
///   outgoing neighbors (v, w).
/// - Perform a bounded BFS from v and from w (depth limited by `max_bubble_len`) to
///   discover meeting nodes m where the two searches converge.
/// - For the best meeting node (smallest combined depth), reconstruct the two
///   u->...->m paths, compute a simple path score (sum of edge lengths) and pick
///   the higher-scoring path to keep. If the higher score is at least
///   `min_support_ratio * lower_score` (e.g. 1.1 to require 10% stronger), then the
///   lower-scoring path is removed (internal nodes removed, excluding u and m).
/// - Removals are RC-aware: the reverse-complement node for each removed node is
///   also removed, and all incoming edges to removed nodes are purged.
///
/// The implementation is deliberately conservative and easy to reason about.

/// Compute reverse-complement oriented node id ("read+" <-> "read-").
fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len()-1];
            format!("{}-", base)
        } else if last == '-' {
            let base = &id[..id.len()-1];
            format!("{}+", base)
        } else {
            id.to_string()
        }
    } else {
        id.to_string()
    }
}

/// Bounded BFS from a start node. Returns maps:
/// - parent: node -> predecessor node (None for the start node)
/// - depth: node -> depth (start has depth 0)
/// - score: node -> summed-edge-length from start to that node (u64)
fn bfs_limited(
    graph: &crate::create_overlap_graph::OverlapGraph,
    start: &str,
    max_depth: usize,
) -> (
    HashMap<String, Option<String>>,
    HashMap<String, usize>,
    HashMap<String, u64>,
) {
    let mut parent: HashMap<String, Option<String>> = HashMap::new();
    let mut depth: HashMap<String, usize> = HashMap::new();
    let mut score: HashMap<String, u64> = HashMap::new();

    let mut q: VecDeque<String> = VecDeque::new();
    parent.insert(start.to_string(), None);
    depth.insert(start.to_string(), 0);
    score.insert(start.to_string(), 0u64);
    q.push_back(start.to_string());

    while let Some(cur) = q.pop_front() {
        let cur_depth = *depth.get(&cur).unwrap_or(&0);
        if cur_depth >= max_depth {
            continue;
        }
        // expand neighbours
        if let Some(node) = graph.nodes.get(&cur) {
            for (tgt, len) in &node.edges {
                // if unseen, record parent/depth/score and enqueue
                if !depth.contains_key(tgt) {
                    parent.insert(tgt.clone(), Some(cur.clone()));
                    depth.insert(tgt.clone(), cur_depth + 1);
                    let parent_score = *score.get(&cur).unwrap_or(&0);
                    score.insert(tgt.clone(), parent_score + (*len as u64));
                    q.push_back(tgt.clone());
                }
            }
        }
    }

    (parent, depth, score)
}

/// Reconstruct path from `start` to `target` using the parent map returned by `bfs_limited`.
/// If `target` is not reachable, returns an empty Vec.
fn reconstruct_path(
    parent: &HashMap<String, Option<String>>,
    start: &str,
    target: &str,
) -> Vec<String> {
    let mut path_rev: Vec<String> = Vec::new();
    let mut cur = target.to_string();
    path_rev.push(cur.clone());
    while let Some(Some(p)) = parent.get(&cur) {
        cur = p.clone();
        path_rev.push(cur.clone());
    }
    // We expect the last item to be start
    if path_rev.last().map(|s| s.as_str()) != Some(start) {
        return Vec::new();
    }
    path_rev.reverse();
    path_rev
}

/// Remove nodes and their reverse-complements from the graph, and purge any edges
/// pointing to removed nodes. Returns number of nodes removed (including RCs).
fn remove_nodes_rc_aware(
    graph: &mut crate::create_overlap_graph::OverlapGraph,
    to_remove: &HashSet<String>,
) -> usize {
    let mut expanded: HashSet<String> = HashSet::new();
    for id in to_remove.iter() {
        expanded.insert(id.clone());
        expanded.insert(rc_node(id));
    }

    // remove nodes
    for id in expanded.iter() {
        graph.nodes.remove(id);
    }

    // purge incoming edges to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges.retain(|(tgt, _)| !expanded.contains(tgt));
    }

    expanded.len()
}

/// Remove simple bubbles in the overlap graph.
///
/// Parameters:
/// - graph: the overlap graph to modify in-place
/// - max_bubble_len: maximum number of edges per side to consider when searching for a converging node
/// - min_support_ratio: minimum ratio required for the stronger path to be considered significantly better
///   e.g. 1.1 means the stronger path must have >= 10% higher score than the weaker one.
///
/// Returns number of oriented nodes removed (counts both nodes and their RC counterparts).
pub fn remove_bubbles(
    graph: &mut crate::create_overlap_graph::OverlapGraph,
    max_bubble_len: usize,
    min_support_ratio: f64,
) -> usize {
    if max_bubble_len == 0 {
        return 0;
    }

    // Snapshot of nodes to iterate safely
    let node_keys: Vec<String> = graph.nodes.keys().cloned().collect();

    let mut total_removed: usize = 0;

    for u in node_keys.iter() {
        // get outgoing neighbors (clone so we don't borrow across mutation)
        let outs = match graph.nodes.get(u) {
            Some(n) => n.edges.iter().map(|(t, l)| (t.clone(), *l)).collect::<Vec<_>>(),
            None => continue,
        };
        if outs.len() < 2 {
            continue;
        }

        // consider every unordered pair of outgoing neighbors
        for i in 0..outs.len() {
            for j in (i + 1)..outs.len() {
                let start_a = &outs[i].0;
                let start_b = &outs[j].0;
                if start_a == start_b {
                    continue;
                }

                let (parent_a, depth_a, score_a) = bfs_limited(graph, start_a, max_bubble_len);
                let (parent_b, depth_b, score_b) = bfs_limited(graph, start_b, max_bubble_len);

                // find meeting nodes
                let reached_a: HashSet<&String> = depth_a.keys().collect();
                let reached_b: HashSet<&String> = depth_b.keys().collect();

                let mut meetings: Vec<(&String, usize)> = Vec::new(); // (node, combined_depth)
                for node in reached_a.intersection(&reached_b) {
                    let d = depth_a.get(*node).unwrap_or(&usize::MAX) + depth_b.get(*node).unwrap_or(&usize::MAX);
                    meetings.push((node.clone(), d));
                }

                if meetings.is_empty() {
                    continue;
                }

                // pick best meeting: minimal combined depth
                meetings.sort_unstable_by_key(|k| k.1);
                let (meet_node, _meet_depth) = meetings[0];

                // reconstruct paths start_a -> meet_node and start_b -> meet_node
                let path_a = reconstruct_path(&parent_a, start_a, meet_node);
                let path_b = reconstruct_path(&parent_b, start_b, meet_node);
                if path_a.is_empty() || path_b.is_empty() {
                    continue;
                }

                // compute path scores (sum of edge lengths) using the BFS score maps
                let score_path_a = *score_a.get(meet_node).unwrap_or(&0u64);
                let score_path_b = *score_b.get(meet_node).unwrap_or(&0u64);

                // If both scores are zero (unexpected), skip
                if score_path_a == 0 && score_path_b == 0 {
                    continue;
                }

                // identify winner and loser
                let (loser_path, loser_score, winner_score) = if score_path_a < score_path_b {
                    (path_a.clone(), score_path_a, score_path_b)
                } else {
                    (path_b.clone(), score_path_b, score_path_a)
                };

                // require winner sufficiently stronger than loser
                if (winner_score as f64) < (min_support_ratio * (loser_score as f64)) {
                    // not a clear bubble by score
                    continue;
                }

                // Nodes to remove: all nodes on loser_path excluding the meeting node
                // Also exclude node `u` (divergence point). Typically the path starts at the neighbor of u
                let mut to_remove: HashSet<String> = HashSet::new();
                for node in loser_path.into_iter() {
                    if node == *meet_node {
                        break;
                    }
                    // Defensive: don't remove u
                    if node == *u {
                        continue;
                    }
                    to_remove.insert(node);
                }

                if to_remove.is_empty() {
                    continue;
                }

                // perform RC-aware removal
                let removed = remove_nodes_rc_aware(graph, &to_remove);
                total_removed += removed;

                // After a modification, the graph changed; to be conservative we stop exploring further
                // pairs from this u, and continue with next u. This avoids complicated borrow/mutation interactions
                // and keeps the algorithm simple and conservative.
                break;
            }
        }
    }

    total_removed
}
