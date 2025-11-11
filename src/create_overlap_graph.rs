/// Overlap graph creation module
/// 1. read filtered PAF alignments
/// 2. classify overlaps into internal matches, containments, proper overlaps
/// 3. build bigraph with nodes for each read in both orientations and directed edges for proper overlaps

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Edge info containing all the metrics we track
#[derive(Clone)]
pub struct EdgeInfo {
    pub target_id: String,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

/// A node in the overlap graph. Earch read is represented by two nodes: "<read_name>+" and "<read_name>-"
/// One for the origininal orientation and one for the reverse complement
/// Each node has directed edges to other nodes with associated edge lengths
pub struct Node {
    node_id: String,
    pub edges: Vec<EdgeInfo>,
}

impl Node {

    /// Create a new node with given id and no edges
    fn new(node_id: String) -> Self {
        Self {
            node_id,
            edges: Vec::new(),
        }
    }

    /// Add a directed edge to a node, if an edge to the target already exists, we ignore (avoid duplicates)
    fn add_edge(&mut self, target_node: &str, edge_len: u32, overlap_len: u32, identity: f64) {
        if self.edges.iter().any(|e| e.target_id == target_node) {
            // multiple edges to the same target
            // currently impossible due to the overlap filtering
            // handling might change in the future to handle this case
            // silently ignore duplicate edges but log for debugging
            eprintln!("Warning: duplicate edge {} -> {} ignored", self.node_id, target_node);
            return;
        }
        self.edges.push(EdgeInfo {target_id: target_node.to_owned(), edge_len, overlap_len, identity,});
    }

    /// Remove a directed edge to the node with target_node id
    pub fn remove_edge(&mut self, target_node: &str) {
        if let Some(pos) = self.edges.iter().position(|e| e.target_id == target_node) {
            self.edges.swap_remove(pos);
        }
    }

    /// Sort edges by length (ascending).
    pub fn sort_edges(&mut self) {
        self.edges.sort_unstable_by_key(|e| e.edge_len);
    }
}

/// Overlap graph containing nodes keyed by their node id
pub struct OverlapGraph {
    pub nodes: HashMap<String, Node>,
}

impl OverlapGraph {

    /// Create a new empty overlap graph
    fn new() -> Self {
        Self {
            nodes: HashMap::new(),
        }
    }

    /// Add a node to the graph if it does not already exist, if it already exists do nothing
    fn add_node(&mut self, node_id: String) {
        self.nodes.entry(node_id.clone()).or_insert_with(|| Node::new(node_id));
    }

    /// Add a directed edge from from_id to to_id with given edge length and metrics
    fn add_edge(&mut self, from_id: &str, to_id: &str, edge_len: u32, overlap_len: u32, identity: f64) {
        // ensure nodes exist
        if !self.nodes.contains_key(from_id) || !self.nodes.contains_key(to_id) {
            panic!("add_edge: nodes must exist before adding edge");
        }
        if let Some(node) = self.nodes.get_mut(from_id) {
            node.add_edge(to_id, edge_len, overlap_len, identity);
        }
    }

    /// Remove a node and all its edges from the graph.
    fn remove_node(&mut self, node_id: &str) {
        self.nodes.remove(node_id);
    }
}

/// Counters returned for diagnostics
#[derive(Default)]
pub struct OverlapGraphStats {
    pub nr_internal_match: u64,
    pub nr_first_contained: u64,
    pub nr_second_contained: u64,
    pub nr_proper_overlaps: u64,
}

/// Parse a PAF line
pub fn parse_paf_line(line: &str) -> Option<(String, u32, i64, i64, char, String, u32, i64, i64, u32, u32, i32)> {
    
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return None;
    }

    // fields indexing:
    // 0: query_name, 1: query_length, 2: query_start, 3: query_end, 4: strand, 5: target_name, 6: target_length, 7: target_start, 8: target_end, 9: num_matching, 10: alignment_block_length, 11: mapq
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

/// Build overlap graph from overlaps
pub fn create_overlap_graph(paf_path: &str, max_overhang: u32, overhang_ratio: f64) -> Result<(OverlapGraph, OverlapGraphStats), std::io::Error> {
    
    // read PAF file
    let reader = BufReader::new(File::open(paf_path)?);

    let mut g = OverlapGraph::new();
    let mut stats = OverlapGraphStats::default();

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

    let (query_name, query_length, query_start, query_end, strand, target_name, target_length, target_start, target_end, num_matching, alignment_block_length, _mapq) = parsed;

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
                let overlap_len = overlap_length as u32;
                let identity = if alignment_block_length == 0 { 0.0 } else { num_matching as f64 / alignment_block_length as f64 };
                g.add_edge(&q_plus, &t_orient, edge1_len, overlap_len, identity);
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
                let overlap_len = overlap_length as u32;
                let identity = if alignment_block_length == 0 { 0.0 } else { num_matching as f64 / alignment_block_length as f64 };
                // direction: t_rc -> q_minus
                g.add_edge(&t_rc, &q_minus, edge2_len, overlap_len, identity);
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
                let overlap_len = overlap_length as u32;
                let identity = if alignment_block_length == 0 { 0.0 } else { num_matching as f64 / alignment_block_length as f64 };
                g.add_edge(&t_orient, &q_plus, edge1_len, overlap_len, identity);
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
                let overlap_len = overlap_length as u32;
                let identity = if alignment_block_length == 0 { 0.0 } else { num_matching as f64 / alignment_block_length as f64 };
                g.add_edge(&q_minus, &t_rc, edge2_len, overlap_len, identity);
            }
        }
    }

    Ok((g, stats))
}