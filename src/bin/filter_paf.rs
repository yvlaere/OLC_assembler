//! A simple tool to filter PAF files based on overlap statistics.
//! Uses a two-pass approach:
//! 1) Collect basic stats per read (number of overlaps, total aligned bases).
//! 2) Filter overlaps based on criteria (e.g., remove self-alignments).

use std::env;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> io::Result<()> {

    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();

    // Check for correct number of arguments
    if args.len() != 3 {
        eprintln!("Usage: {} <input.paf> <output.paf>", args[0]);
        std::process::exit(1);
    }

    // Input and output file paths
    let input_path = &args[1];
    let output_path = &args[2];

    // Open input file for reading
    let reader = BufReader::new(File::open(input_path)?);
    let mut output_file = File::create(output_path)?;

    // Declaring variables
    // Mapping from sequence names to IDs, the ID is the index in the stats vectors (so they are integers starting from 0, incrementing by 1)
    let mut name2id: HashMap<String, usize> = HashMap::new();
    let mut next_id: usize = 0;
    // Stats vectors, indexed by the sequence ID
    let mut overlap_count: Vec<u32> = Vec::new();
    let mut sum_aligned_bases: Vec<u64> = Vec::new();
    let mut read_length: Vec<u32> = Vec::new();
    let min_overlap_lenth = 500;

    // First pass: gather statistics
    for line_res in reader.lines() {

        // Skip comments and empty lines and split fields
        let line = line_res?;
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 { continue; }

        // Parse sequence names
        let query_name = fields[0];
        let target_name = fields[5];

        // Skip self-alignments
        if query_name == target_name { continue; }

        // Parse remaining fields
        let query_length: usize = fields[1].parse().unwrap_or(0);
        let query_start: usize = fields[2].parse().unwrap_or(0);
        let query_end: usize = fields[3].parse().unwrap_or(0);
        let strand = fields[4];
        let target_length: usize = fields[6].parse().unwrap_or(0);
        let target_start: usize = fields[7].parse().unwrap_or(0);
        let target_end: usize = fields[8].parse().unwrap_or(0);
        let num_matching: usize = fields[9].parse().unwrap_or(0);
        let alignment_block_length: usize = fields[10].parse().unwrap_or(0);
        let mapping_quality: usize = fields[11].parse().unwrap_or(0);

        // Get the query and target IDs from their names, creating a new ID if it they didn’t exist yet
        // and initializing their stat
        let query_id = *name2id.entry(query_name.clone()).or_insert_with(|| {
            // Create a new ID
            let id = next_id;
            next_id += 1;
            overlap_count.push(0);
            sum_aligned_bases.push(0);
            read_length.push(query_length);
            // Retrun the new ID
            id
        });
        let target_id = *name2id.entry(target_name.clone()).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            overlap_count.push(0);
            sum_aligned_bases.push(0);
            read_length.push(target_length);
            id
        });

        // Only count overlaps passing some minimal filters:
        if alignment_block_length < min_overlap_lenth { continue; }

        overlap_count[query_id] += 1;
        overlap_count[target_id] += 1;
        sum_aligned_bases[query_id] += alignment_block_length as u64;
        sum_aligned_bases[target_id] += alignment_block_length as u64;
    }

    Ok(())
}