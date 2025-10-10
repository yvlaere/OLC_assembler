use std::env; // For command-line argument parsing
use std::fs::File; // For file handling
use std::io::{self, BufRead, BufReader, Write}; // For I/O operations

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
    let input_file = File::open(input_path)?;
    let reader = BufReader::new(input_file);
    let mut output_file = File::create(output_path)?;

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        // Skip malformed lines
        if fields.len() < 12 {
            continue;
        }

        // Parse necessary fields
        let query_name = fields[0];
        let target_name = fields[5];

        // Skip self-alignments
        if query_name == target_name {
            continue; // Skip self-alignments
        }

        let query_start: usize = fields[2].parse().unwrap_or(0);
        let query_end: usize = fields[3].parse().unwrap_or(0);
        let target_start: usize = fields[7].parse().unwrap_or(0);
        let target_end: usize = fields[8].parse().unwrap_or(0);

        //if (query_end > query_start) && (target_end > target_start) {
        //    writeln!(output_file, "{}", line)?;
        //}
    }

    Ok(())
}