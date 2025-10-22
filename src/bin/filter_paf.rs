//! A tool to filter PAF files based on overlap statistics.
//! Uses a two-pass approach:
//! 1) Collect basic stats per read (number of overlaps, aligned bases)
//! 2) Filter overlaps based on quality criteria

use std::env;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

struct PafRecord {
    query_name: String,
    query_length: u32,
    target_name: String,
    target_length: u32,
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
            target_name: fields[5].to_string(),
            target_length: fields[6].parse().ok()?,
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

fn main() -> io::Result<()> {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.paf> <output.paf>", args[0]);
        std::process::exit(1);
    }

    // Configuration
    let min_overlap_length: u32 = 2000;
    let min_overlap_count: u32 = 3;
    let min_percent_identity: f32 = 85.0;

    // Setup data structures
    let mut name2id: HashMap<String, usize> = HashMap::new();
    let mut overlap_counts: Vec<u32> = Vec::new();
    let mut ids2overlap: HashMap<(usize, usize), u32> = HashMap::new();
    let mut next_id: usize = 0;

    // First pass: collect statistics
    let reader = BufReader::new(File::open(&args[1])?);
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        if let Some(record) = PafRecord::from_line(&line) {
            if record.is_self_alignment() { continue; }
            
            // Skip if it doesn't meet basic criteria
            if record.alignment_block_length < min_overlap_length || 
               record.percent_identity() < min_percent_identity {
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

            // Update statistics
            overlap_counts[query_id] += 1;
            overlap_counts[target_id] += 1;

            let key = if query_id < target_id {
                (query_id, target_id)
            } else {
                (target_id, query_id)
            };

            ids2overlap
                .entry(key)
                .and_modify(|e| *e = (*e).max(record.alignment_block_length))
                .or_insert(record.alignment_block_length);
        }
    }

    // Second pass: write filtered records
    let reader = BufReader::new(File::open(&args[1])?);
    let mut writer = File::create(&args[2])?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }

        if let Some(record) = PafRecord::from_line(&line) {
            if record.is_self_alignment() { continue; }

            if let (Some(&query_id), Some(&target_id)) = (name2id.get(&record.query_name), name2id.get(&record.target_name)) {
                if overlap_counts[query_id] >= min_overlap_count && 
                   overlap_counts[target_id] >= min_overlap_count &&
                   record.alignment_block_length >= min_overlap_length &&
                   record.percent_identity() >= min_percent_identity {
                    
                    let key = if query_id < target_id {
                        (query_id, target_id)
                    } else {
                        (target_id, query_id)
                    };

                    if let Some(&best_len) = ids2overlap.get(&key) {
                        if record.alignment_block_length == best_len {
                            writeln!(writer, "{}", line)?;
                        }
                    }
                }
            }
        }
    }

    Ok(())
}