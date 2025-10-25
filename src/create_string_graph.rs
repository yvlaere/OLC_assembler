//! A tool to create a string graph from PAF overlaps.
//! based on the the miniasm pseudo-code.

use std::env;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// A node in the assembly graph. Earch read is represented by two nodes: "<read_name>+" and "<read_name>-". 
/// One for the origininal orientation and one for the reverse complement.
/// Each node has directed edges to other nodes with associated edge lengths.
struct Node {
    node_id: String,
    // edges: vector of (target_node_id, edge_length)
    edges: Vec<(String, u32)>,
}

impl Node {
    fn new(node_id: String) -> Self {
        Self {
            node_id,
            edges: Vec::new(),
        }
    }

    /// Add a directed edge to this node. If an edge to the target already exists, we ignore (avoid duplicates).
    fn add_edge(&mut self, target_node: &str, edge_len: u32) {
        if self.edges.iter().any(|(t, _)| t == target_node) {
            // multiple edges to the same target
            // currently impossible due to the overlap filtering
            // handling might change in the future to handle this case
            return;
        }
        self.edges.push((target_node.to_owned(), edge_len));
    }

    /// Remove a directed edge to the node with target_node id.
    fn remove_edge(&mut self, target_node: &str) {
        if let Some(pos) = self.edges.iter().position(|(t, _)| t == target_node) {
            self.edges.swap_remove(pos);
        }
    }

    /// Sort edges by length (ascending).
    fn sort_edges(&mut self) {
        self.edges.sort_unstable_by_key(|(_, len)| *len);
    }
}

/// Assembly graph containing nodes keyed by "<read_name><+/->"
struct AssemblyGraph {
    nodes: HashMap<String, Node>,
}

impl AssemblyGraph {
    fn new() -> Self {
        Self {
            nodes: HashMap::new(),
        }
    }

    /// Add a node to the graph if it does not already exist, if it already exists do nothing.
    fn add_node(&mut self, node_id: String) {
        self.nodes.entry(node_id.clone()).or_insert_with(|| Node::new(node_id));
    }

    /// Add a directed edge from from_id to to_id with given edge length.
    fn add_edge(&mut self, from_id: &str, to_id: &str, edge_len: u32) {
        // ensure nodes exist
        if !self.nodes.contains_key(from_id) || !self.nodes.contains_key(to_id) {
            panic!("add_edge: nodes must exist before adding edge");
        }
        if let Some(node) = self.nodes.get_mut(from_id) {
            node.add_edge(to_id, edge_len);
        }
    }

    /// Remove a node and all its edges from the graph.
    fn remove_node(&mut self, node_id: &str) {
        self.nodes.remove(node_id);
    }

    // Sort edges in all nodes by length.
    //fn sort_all_edges(&mut self) {
    //    for node in self.nodes.values_mut() {
    //        node.sort_edges();
    //    }
    //}
}

/// Counters returned for diagnostics
#[derive(Default)]
struct GraphStats {
    nr_internal_match: u64,
    nr_first_contained: u64,
    nr_second_contained: u64,
    nr_proper_overlaps: u64,
}

/// Parse a PAF line
fn parse_paf_line(line: &str) -> Option<(String, u32, i64, i64, char, String, u32, i64, i64, u32, u32, i32)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return None;
    }

    // fields indexing:
    // 0: query_name, 1: query_length, 2: query_start, 3: query_end, 4: strand
    // 5: target_name, 6: target_length, 7: target_start, 8: target_end, 9: num_matching, 10: alignment_block_length, 11: mapq
    let query_name = fields[0].to_string();
    let query_length = fields[1].parse::<u32>().ok()?;
    let query_start = fields[2].parse::<i64>().ok()?;
    let query_end = fields[3].parse::<i64>().ok()?;
    let strand = fields[4].chars().next().unwrap_or('+');
    let target_name = fields[5].to_string();
    let target_length = fields[6].parse::<u32>().ok()?;
    let target_start = fields[7].parse::<i64>().ok()?;
    let target_end = fields[8].parse::<i64>().ok()?;
    let num_matching = fields[9].parse::<u32>().ok()?;
    let alignment_block_length = fields[10].parse::<u32>().ok()?;
    let mapq = fields[11].parse::<i32>().unwrap_or(0);
    Some((query_name, query_length, query_start, query_end, strand, target_name, target_length, target_start, target_end, num_matching, alignment_block_length, mapq))
}

/// Build string graph from overlaps
fn create_string_graph(paf_path: &str, max_overhang: u32, overhang_ratio: f64) -> Result<(AssemblyGraph, GraphStats), std::io::Error> {
    
    // read PAF file
    let reader = BufReader::new(File::open(paf_path)?);

    let mut g = AssemblyGraph::new();
    let mut stats = GraphStats::default();

    for (line_nr, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        let parsed = match parse_paf_line(&line) {
            Some(p) => p,
            None => {
                // skip malformed lines
                eprintln!("Warning: skipping malformed PAF line {}", line_nr + 1);
                continue;
            }
        };

        let (query_name, query_length, query_start, query_end, strand, target_name, target_length, target_start, target_end, _num_matching, _alignment_block_length, _mapq) = parsed;

        // use i64 for arithmetic because we may subtract and want to allow signed intermediate
        let b1 = query_start;
        let e1 = query_end;
        let l1 = query_length as i64;

        // define overlap beginning and end based on orientation
        // using naming convention corresponding with miniasm paper
        let (b2, e2, l2) = if strand == '+' {
            (target_start, target_end, target_length as i64)
        } else {
            // reverse complement coordinates on the target
            (target_length as i64 - target_end, target_length as i64 - target_start, target_length as i64)
        };

        // overhang is the part next to the overlap where the reads don't align, but should in case of perfect overlap
        // overhang: min(b1,b2) + min(l1 - e1, l2 - e2)
        let overhang_left = std::cmp::min(b1, b2);
        let overhang_right = std::cmp::min(l1 - e1, l2 - e2);
        let overhang = overhang_left + overhang_right;

        // longest overlap length (max aligned spans on either read)
        let overlap_length1 = e1.saturating_sub(b1);
        let overlap_length2 = e2.saturating_sub(b2);
        let overlap_length = std::cmp::max(overlap_length1, overlap_length2) as f64;

        // decide overhang threshold: max_overhang (u32) or maplen * overhang_ratio
        let overhang_threshold = (overlap_length * overhang_ratio).ceil() as i64;
        let max_overhang_i64 = max_overhang as i64;
        let allowed_overhang = std::cmp::min(max_overhang_i64, overhang_threshold);

        // classification of the overlap
        // overlaps are a subset of alignments where (in theory) two read edges, one from each read, are part of the alignment
        if overhang > allowed_overhang {
            // internal match
            stats.nr_internal_match += 1;
            continue;
        }

        // conditions for containment:
        // first contained in second:
        // (b1 <= b2) and ((l1 - e1) <= (l2 - e2))
        let first_contained = (b1 <= b2) && ((l1 - e1) <= (l2 - e2));
        // second contained in first:
        let second_contained = (b1 >= b2) && ((l1 - e1) >= (l2 - e2));

        if first_contained {
            stats.nr_first_contained += 1;
            continue;
        } else if second_contained {
            stats.nr_second_contained += 1;
            continue;
        }

        // at this point it's a proper overlap between reads
        // decide orientation and edge lengths
        if b1 > b2 {
            // first read to second read overlap (query -> target)
            stats.nr_proper_overlaps += 1;

            // add nodes & edges (forward orientation)
            let q_plus = format!("{}+", query_name);
            let t_orient = format!("{}{}", target_name, strand);
            g.add_node(q_plus.clone());
            g.add_node(t_orient.clone());

            // edge length = b1 - b2 (non-overlapping prefix length)
            let edge1_len_i64 = b1 - b2;
            if edge1_len_i64 >= 0 {
                let edge1_len = edge1_len_i64 as u32;
                g.add_edge(&q_plus, &t_orient, edge1_len);
            }

            // reverse complement counterpart:
            let q_minus = format!("{}-", query_name);
            let rc_strand = if strand == '+' { '-' } else { '+' };
            let t_rc = format!("{}{}", target_name, rc_strand);
            g.add_node(q_minus.clone());
            g.add_node(t_rc.clone());

            // edge2_len = (l2 - e2) - (l1 - e1)
            let edge2_len_i64 = (l2 - e2) - (l1 - e1);
            if edge2_len_i64 >= 0 {
                let edge2_len = edge2_len_i64 as u32;
                // direction: t_rc -> q_minus
                g.add_edge(&t_rc, &q_minus, edge2_len);
            }
        } else {
            // second to first overlap (target -> query)
            stats.nr_proper_overlaps += 1;

            let q_plus = format!("{}+", query_name);
            let t_orient = format!("{}{}", target_name, strand);
            g.add_node(q_plus.clone());
            g.add_node(t_orient.clone());

            let edge1_len_i64 = b2 - b1;
            if edge1_len_i64 >= 0 {
                let edge1_len = edge1_len_i64 as u32;
                // direction t -> q
                g.add_edge(&t_orient, &q_plus, edge1_len);
            }

            // reverse complement counterpart:
            let q_minus = format!("{}-", query_name);
            let rc_strand = if strand == '+' { '-' } else { '+' };
            let t_rc = format!("{}{}", target_name, rc_strand);
            g.add_node(q_minus.clone());
            g.add_node(t_rc.clone());

            let edge2_len_i64 = (l1 - e1) - (l2 - e2);
            if edge2_len_i64 >= 0 {
                let edge2_len = edge2_len_i64 as u32;
                // direction q_minus -> t_rc
                g.add_edge(&q_minus, &t_rc, edge2_len);
            }
        }
    }

    Ok((g, stats))
}