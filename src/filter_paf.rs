/// PAF filtering module
/// Uses a two-pass approach:
/// 1) Collect basic stats per read (number of overlaps, aligned bases)
/// 2) Filter overlaps based on quality criteria

use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Struct to hold a PAF record
struct PafRecord {
    query_name: String,
    query_length: u32,
    query_start: i64,
    query_end: i64,
    strand: char,
    target_name: String,
    target_length: u32,
    target_start: i64,
    target_end: i64,
    num_matching: u32,
    alignment_block_length: u32,
}

// Implement methods for PafRecord (kind of like class methods)
impl PafRecord {
    fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 { return None; }
        
        Some(Self {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().ok()?,
            query_start: fields[2].parse::<i64>().ok()?,
            query_end: fields[3].parse::<i64>().ok()?,
            strand: fields[4].chars().next().unwrap_or('+'),
            target_name: fields[5].to_string(),
            target_length: fields[6].parse::<u32>().ok()?,
            target_start: fields[7].parse::<i64>().ok()?,
            target_end: fields[8].parse::<i64>().ok()?,
            num_matching: fields[9].parse().ok()?,
            alignment_block_length: fields[10].parse().ok()?,
        })
    }

    fn is_self_alignment(&self) -> bool {
        self.query_name == self.target_name
    }

    fn percent_identity(&self) -> f32 {
        (self.num_matching as f32 / self.alignment_block_length as f32) * 100.0
    }
}

/// Filter PAF file based on overlap quality criteria
pub fn filter_paf(paf_in: &str, paf_out: &str, min_overlap_length: &u32, min_overlap_count: &u32, min_percent_identity: &f32, max_overhang: &u32, overhang_ratio: &f64) -> io::Result<()> {

    // Setup data structures
    let mut name2id: HashMap<String, usize> = HashMap::new();
    let mut overlap_counts: Vec<u32> = Vec::new();
    let mut ids2overlap: HashMap<(usize, usize), u32> = HashMap::new();
    let mut next_id: usize = 0;

    // NEW: track reads that are contained in another read (we will exclude them)
    let mut contained_reads: HashSet<usize> = HashSet::new();

    // First pass: collect statistics
    let reader = BufReader::new(File::open(paf_in)?);
    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') || line.trim().is_empty() { continue; }

        if let Some(record) = PafRecord::from_line(&line) {

            // skip self-alignments
            if record.is_self_alignment() { continue; }
            
            // skip short alignment
            if record.alignment_block_length < *min_overlap_length {
                continue;
            }

            // skip low identity
            if record.percent_identity() < *min_percent_identity {
                continue;
            }

            // Get or create IDs for query and target
            let query_id = match name2id.get(&record.query_name) {
                Some(&id) => id,
                None => {
                    let id = next_id;
                    next_id += 1;
                    overlap_counts.push(0);
                    name2id.insert(record.query_name.clone(), id);
                    id
                }
            };

            let target_id = match name2id.get(&record.target_name) {
                Some(&id) => id,
                None => {
                    let id = next_id;
                    next_id += 1;
                    overlap_counts.push(0);
                    name2id.insert(record.target_name.clone(), id);
                    id
                }
            };

            // skip contained
            let b1 = record.query_start;
            let e1 = record.query_end;
            let l1 = record.query_length as i64;

            // define overlap beginning and end based on orientation
            // using naming convention corresponding with miniasm paper
            let (b2, e2, l2) = if record.strand == '+' {
                (record.target_start, record.target_end, record.target_length as i64)
            } else {
                // reverse complement coordinates on the target
                (record.target_length as i64 - record.target_end, record.target_length as i64 - record.target_start, record.target_length as i64)
            };

            // first contained in second:
            // (b1 <= b2) and ((l1 - e1) <= (l2 - e2))
            let first_contained = (b1 <= b2) && ((l1 - e1) <= (l2 - e2));
            // second contained in first:
            let second_contained = (b1 >= b2) && ((l1 - e1) >= (l2 - e2));
            if first_contained {
                contained_reads.insert(query_id);
                continue;
            } else if second_contained {
                contained_reads.insert(target_id);
                continue;
            }

            // Update statistics
            overlap_counts[query_id] += 1;
            overlap_counts[target_id] += 1;

            let key = if query_id < target_id {
                (query_id, target_id)
            } else {
                (target_id, query_id)
            };

            ids2overlap.entry(key).and_modify(|e| *e = (*e).max(record.alignment_block_length)).or_insert(record.alignment_block_length);
        }
    }

    // Second pass: write filtered records
    let reader = BufReader::new(File::open(paf_in)?);
    let mut writer = File::create(paf_out)?;

    let mut self_alignments_skipped = 0;
    let mut contained_skipped = 0;
    let mut overlap_count_skipped = 0;
    let mut alignment_length_skipped = 0;
    let mut percent_identity_skipped = 0;
    let mut best_overlap_skipped = 0;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        if let Some(record) = PafRecord::from_line(&line) {
            if record.is_self_alignment() { self_alignments_skipped +=1; continue; }

            if let (Some(&query_id), Some(&target_id)) = (name2id.get(&record.query_name), name2id.get(&record.target_name)) {
                
                if record.alignment_block_length < *min_overlap_length {
                    alignment_length_skipped += 1;
                    continue;
                }

                if record.percent_identity() < *min_percent_identity {
                    percent_identity_skipped += 1;
                    continue;
                }
                
                // NEW: skip any record involving a contained read
                if contained_reads.contains(&query_id) || contained_reads.contains(&target_id) {
                    contained_skipped += 1;
                    continue;
                }
                
                if overlap_counts[query_id] < *min_overlap_count || overlap_counts[target_id] < *min_overlap_count {
                    overlap_count_skipped += 1;
                    continue;
                }

                
                
                let key = if query_id < target_id {
                    (query_id, target_id)
                } 
                else {
                    (target_id, query_id)
                };

                if let Some(&best_len) = ids2overlap.get(&key) {
                    if record.alignment_block_length == best_len {
                        writeln!(writer, "{}", line)?;
                    }
                    else {
                        best_overlap_skipped += 1;
                        // not the best overlap, skip
                    }
                }
                
            }
        }
    }

    println!("PAF filtering summary:");
    println!("  Self-alignments skipped: {}", self_alignments_skipped);
    println!("  Contained reads skipped: {}", contained_skipped);
    println!("  Overlap count filter skipped: {}", overlap_count_skipped);
    println!("  Alignment length filter skipped: {}", alignment_length_skipped);
    println!("  Percent identity filter skipped: {}", percent_identity_skipped);
    println!("  Best overlap filter skipped: {}", best_overlap_skipped);

    Ok(())
}