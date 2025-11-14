/// alignment filtering module
/// Uses a two-pass approach:
/// 1) Read all alignments, store them if they pass basic filters (self-alignment, overlap length, identity), only store the longest alignment per read pair
/// 2) Calculate coverage statistics per reads
/// 3) Classify alignments into internal matches, contained reads, proper overlaps
/// 4) Filter the proper overlaps

use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Struct to hold a read
struct Read {
    id: usize,
    name: String,
    length: u32,
    per_base_coverage: Vec<u32>,
}

/// Struct to hold an alignment
#[derive(Clone)]
struct Alignment {
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
    mapq: u8,
}

// Struct to hold overlap information
pub struct Overlap {
    pub query_name: String,
    pub rc_query_name: String,
    pub target_name: String,
    pub rc_target_name: String,
    pub edge_len: u32,
    pub rc_edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

/// Enum for alignment classification
enum AlignmentType {
    InternalMatch,
    FirstContained,
    SecondContained,
    ProperOverlap,
}

// Implement methods for Alignment (kind of like class methods)
impl Alignment {
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
            mapq: fields[11].parse().ok()?,
        })
    }

    fn is_self_alignment(&self) -> bool {
        self.query_name == self.target_name
    }

    fn percent_identity(&self) -> f32 {
        (self.num_matching as f32 / self.alignment_block_length as f32) * 100.0
    }
}

/// Classify alignment and update contained reads set
fn classify_alignment(r: &Alignment, query_id: usize, target_id: usize, overlaps: &mut HashMap<(usize, usize), Overlap>, max_overhang: u32, overhang_ratio: f64) -> AlignmentType {
    
    // use i64 for arithmetic because we may subtract and want to allow signed intermediate
    let b1 = r.query_start;
    let e1 = r.query_end;
    let l1 = r.query_length as i64;

    // define overlap beginning and end based on orientation
    // using naming convention corresponding with miniasm paper
    let (b2, e2, l2) = if r.strand == '+' {
        (r.target_start, r.target_end, r.target_length as i64)
    } else {
        // reverse complement coordinates on the target
        (r.target_length as i64 - r.target_end, r.target_length as i64 - r.target_start, r.target_length as i64)
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
        return AlignmentType::InternalMatch;
    }

    // conditions for containment:
    // first contained in second:
    // (b1 <= b2) and ((l1 - e1) <= (l2 - e2))
    let first_contained = (b1 <= b2) && ((l1 - e1) <= (l2 - e2));
    // second contained in first:
    let second_contained = (b1 >= b2) && ((l1 - e1) >= (l2 - e2));

    if first_contained {
        return AlignmentType::FirstContained;
    } else if second_contained {
        return AlignmentType::SecondContained;
    }

    // at this point it's a proper overlap between reads
    // decide orientation and edge lengths
    if b1 > b2 {
        // first read to second read overlap (query -> target)

        // get nodes & edges (forward orientation)
        let q_plus = format!("{}+", r.query_name);
        let t_orient = format!("{}{}", r.target_name, r.strand);

        // edge length = b1 - b2 (non-overlapping prefix length)
        let edge1_len_i64 = b1 - b2;
        if edge1_len_i64 < 0 {
            println!("Warning: negative edge length encountered in overlap classification.");
        }
        let edge1_len = edge1_len_i64 as u32;

        // reverse complement counterpart:
        // direction: t_rc -> q_minus
        let q_minus = format!("{}-", r.query_name);
        let rc_strand = if r.strand == '+' { '-' } else { '+' };
        let t_rc = format!("{}{}", r.target_name, rc_strand);

        // edge length = (l2 - e2) - (l1 - e1)
        let edge2_len_i64 = (l2 - e2) - (l1 - e1);
        if edge2_len_i64 < 0 {
            println!("Warning: negative edge length encountered in overlap classification.");
        }
        let edge2_len = edge2_len_i64 as u32;

        // shared stats
        let overlap_len = r.alignment_block_length as u32;
        let identity = r.percent_identity() as f64;

        // create overlap object and store
        let ov = Overlap {
            query_name: q_plus,
            rc_query_name: q_minus,
            target_name: t_orient,
            rc_target_name: t_rc,
            edge_len: edge1_len as u32,
            rc_edge_len: edge2_len as u32,
            overlap_len: overlap_length as u32,
            identity,
        };
        overlaps.insert((query_id, target_id), ov);

    } 
    else {

        // second to first overlap (target -> query)
        // direction t -> q
        let q_plus = format!("{}+", r.query_name);
        let t_orient = format!("{}{}", r.target_name, r.strand);

        let edge1_len_i64 = b2 - b1;
        if edge1_len_i64 < 0 {
            println!("Warning: negative edge length encountered in overlap classification.");
        }
        let edge1_len = edge1_len_i64 as u32;

        // reverse complement counterpart:
        // direction q_minus -> t_rc
        let q_minus = format!("{}-", r.query_name);
        let rc_strand = if r.strand == '+' { '-' } else { '+' };
        let t_rc = format!("{}{}", r.target_name, rc_strand);

        let edge2_len_i64 = (l1 - e1) - (l2 - e2);
        if edge2_len_i64 < 0 {
            println!("Warning: negative edge length encountered in overlap classification.");
        }
        let edge2_len = edge2_len_i64 as u32;

        // shared stats
        let overlap_len = r.alignment_block_length as u32;
        let identity = r.percent_identity() as f64;

        // create overlap object and store
        let ov = Overlap {
            query_name: q_plus,
            rc_query_name: q_minus,
            target_name: t_orient,
            rc_target_name: t_rc,
            edge_len: edge1_len as u32,
            rc_edge_len: edge2_len as u32,
            overlap_len: overlap_length as u32,
            identity,
        };
        overlaps.insert((query_id, target_id), ov);
    }
    AlignmentType::ProperOverlap
}

/// Filter PAF file based on overlap quality criteria
pub fn filter_paf(paf_in: &str, paf_out: &str, min_overlap_length: &u32, min_overlap_count: &u32, min_percent_identity: &f32, max_overhang: &u32, overhang_ratio: &f64) -> Result<HashMap<(usize, usize), Overlap>, io::Error> {

    // Setup data structures
    // read name to read id mapping
    let mut read_name2read_id: HashMap<String, usize> = HashMap::new();
    // vector of read objects, indexed by read id
    let mut reads: Vec<Read> = Vec::new();
    // maps tuples of read ids to alignment objects
    let mut alignments: HashMap<(usize, usize), Alignment> = HashMap::new();
    // vector to track existing alignments per read id
    // useful for querying alignments
    let mut alignment_ids_per_read: Vec<HashSet<usize>> = Vec::new();
    // initialize read id
    let mut next_id: usize = 0;
    // keep track of contained reads
    let mut contained_reads: HashSet<usize> = HashSet::new();
    // initialize overlap storage
    let mut overlaps: HashMap<(usize, usize), Overlap> = HashMap::new();

    let mut self_alignments_skipped: usize = 0;
    let mut alignment_length_skipped: usize = 0;
    let mut percent_identity_skipped: usize = 0;

    // read the alignments from the PAF file
    let reader = BufReader::new(File::open(paf_in)?);
    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') || line.trim().is_empty() { continue; }

        if let Some(record) = Alignment::from_line(&line) {

            // skip if it's a self-alignment, short alignment, or low identity
            //if record.is_self_alignment() || record.alignment_block_length < *min_overlap_length || record.percent_identity() < *min_percent_identity { continue; }
            if record.is_self_alignment() { self_alignments_skipped +=1; continue; }
            
            if record.alignment_block_length < *min_overlap_length {
                alignment_length_skipped += 1;
                continue;
            }
            if record.percent_identity() < *min_percent_identity {
                percent_identity_skipped += 1;
                continue;
            }
            // Get or create read ids for query and target
            let query_id = match read_name2read_id.get(&record.query_name) {
                
                // we have seen this read before, get its id
                Some(&id) => id,

                // new read, assign new id and create read object
                None => {
                    let id = next_id;
                    next_id += 1;
                    read_name2read_id.insert(record.query_name.clone(), id);
                    alignment_ids_per_read.push(HashSet::new());
                    // create new read object
                    reads.push(Read {id, name: record.query_name.clone(), length: record.query_length, per_base_coverage: vec![0; record.query_length as usize],
                    });
                    id
                }
            };

            let target_id = match read_name2read_id.get(&record.target_name) {
                Some(&id) => id,
                None => {
                    let id = next_id;
                    next_id += 1;
                    read_name2read_id.insert(record.target_name.clone(), id);
                    alignment_ids_per_read.push(HashSet::new());
                    // create new read object
                    reads.push(Read {id, name: record.target_name.clone(), length: record.target_length, per_base_coverage: vec![0; record.target_length as usize],});
                    id
                }
            };

            // extract needed info from the record
            let (qstart, qend) = (record.query_start as usize, record.query_end as usize);
            let (tstart, tend) = (record.target_start as usize, record.target_end as usize);

            // store alignment record
            // if multiple alignments exist between the same read pair, keep the longest one

            if alignment_ids_per_read[query_id].contains(&target_id) {
                
                // an alignment between these reads already exists, it may be stored under
                // (query_id, target_id) or (target_id, query_id)
                if let Some(existing) = alignments.get(&(query_id, target_id)) {
                    if record.alignment_block_length > existing.alignment_block_length {
                        // replace the existing directed entry
                        alignments.insert((query_id, target_id), record.clone());
                    }
                } 
                else if let Some(existing) = alignments.get(&(target_id, query_id)) {
                    if record.alignment_block_length > existing.alignment_block_length {
                        // remove the old reversed entry and store the new (keeps orientation of current record)
                        alignments.remove(&(target_id, query_id));
                        alignments.insert((query_id, target_id), record.clone());
                    }
                } else {
                    println!("Warning: alignment existence inconsistency detected.");
                }
            }
            // we don't have an alignment between these reads yet
            else {
                alignment_ids_per_read[query_id].insert(target_id);
                alignment_ids_per_read[target_id].insert(query_id);
                alignments.insert((query_id, target_id), record);
            }

            // update read coverage statistics
            reads[query_id].per_base_coverage[qstart as usize .. qend as usize].iter_mut().for_each(|c| *c += 1);
            reads[target_id].per_base_coverage[tstart as usize .. tend as usize].iter_mut().for_each(|c| *c += 1);
        }
    }

    println!("Total reads processed: {}", reads.len());
    println!("Total self-alignments skipped: {}", self_alignments_skipped);
    println!("Total alignments skipped due to length filter: {}", alignment_length_skipped);
    println!("Total alignments skipped due to percent identity filter: {}", percent_identity_skipped);

    // count amount of alignments
    println!("Total alignments read: {}", alignments.len());

    // all alignments have been read, now classify them and update contained reads set
    for ((query_id, target_id), alignment) in &alignments {
        let alignment_type = match classify_alignment(alignment, *query_id, *target_id, &mut overlaps, *max_overhang, *overhang_ratio) {
            AlignmentType::InternalMatch => {
                continue; // skip internal matches
            }
            AlignmentType::FirstContained => {
                // mark first read as contained
                contained_reads.insert(*query_id);
                continue;
            }
            AlignmentType::SecondContained => {
                // mark second read as contained
                contained_reads.insert(*target_id);
                continue;
            }
            AlignmentType::ProperOverlap => {
                continue;
            }
        };
    }

    println!("Total overlaps after filtering: {}", overlaps.len());

    // filter contained reads from overlaps
    overlaps.retain(|(q_id, t_id), _| !contained_reads.contains(q_id) && !contained_reads.contains(t_id));

    println!("Total overlaps after removing contained reads: {}", overlaps.len());

    // filter low coverage reads
    let threshold = *min_overlap_count as u32;
    let low_coverage_reads: Vec<_> = reads.iter().enumerate().filter(|(_, r)| r.per_base_coverage.iter().copied().max().unwrap_or(0) < threshold).map(|(id, _)| id).collect();
    overlaps.retain(|(q_id, t_id), _| !low_coverage_reads.contains(q_id) && !low_coverage_reads.contains(t_id));
    println!("Total overlaps after removing low coverage reads: {}", overlaps.len());


    fn per_read_stats(per_base: &[u32]) -> (f64, f64, f64) {
        let len = per_base.len() as f64;

        if len == 0.0 {
            return (0.0, 0.0, 0.0);
        }

        let total: f64 = per_base.iter().map(|&x| x as f64).sum();
        let mean = total / len;

        let cov3 = per_base.iter().filter(|&&x| x >= 3).count() as f64;
        let pct_cov3 = (cov3 / len) * 100.0;

        (mean, pct_cov3, cov3)
    }

    fn global_coverage_stats(reads: &[Read]) {
        let mut global_total_cov: f64 = 0.0;
        let mut global_total_bases: f64 = 0.0;
        let mut global_cov3_bases: f64 = 0.0;

        let mut per_read_means: Vec<f64> = Vec::new();
        let mut per_read_pct3: Vec<f64> = Vec::new();
        let mut per_read_cov3_counts: Vec<f64> = Vec::new();

        let total_nr_reads = reads.len();

        for read in reads {
            let cov = &read.per_base_coverage;
            if cov.is_empty() { continue; }

            // per-read stats
            let (mean, pct3, cov3) = per_read_stats(cov);

            per_read_means.push(mean);
            per_read_pct3.push(pct3);
            per_read_cov3_counts.push(cov3);

            // accumulate global totals
            global_total_bases += cov.len() as f64;
            global_total_cov += cov.iter().map(|&x| x as f64).sum::<f64>();
            global_cov3_bases += cov.iter().filter(|&&x| x >= 3).count() as f64;
        }

        let global_mean = global_total_cov / global_total_bases;
        let global_pct3 = (global_cov3_bases / global_total_bases) * 100.0;
        let mean_cov3 = global_cov3_bases / per_read_pct3.len().max(1) as f64;

        let mean_of_means =
            per_read_means.iter().sum::<f64>() / per_read_means.len().max(1) as f64;

        let mean_of_pct3 =
            per_read_pct3.iter().sum::<f64>() / per_read_pct3.len().max(1) as f64;

        let pct3_total =
            per_read_pct3.len().max(1) as f64 / total_nr_reads.max(1) as f64;

        println!("=== GLOBAL COVERAGE STATISTICS ===");
        println!("Total bases across all reads: {}", global_total_bases as u64);
        println!("Global mean coverage: {:.3}", global_mean);
        println!("Global % bases with cov ≥3: {:.3}%", global_pct3);
        println!("Mean of per-read mean coverages: {:.3}", mean_of_means);
        println!("Mean of per-read pct ≥3: {:.3}%", mean_of_pct3);
        println!("Fraction of reads with any bases cov ≥3: {:.3}%", pct3_total * 100.0);
        println!("Mean # bases with cov ≥3 per read: {:.3}", mean_cov3);
    }

    global_coverage_stats(&reads);


    // write filtered overlaps to output PAF file
    let mut writer = io::BufWriter::new(File::create(paf_out)?);
    for ((_q_id, _t_id), ov) in &overlaps {
        // write the corresponding Alignment record to a new paf file
        // start by getting the corresponding Alignment object
        let alignment = alignments.get(&(*_q_id, *_t_id)).unwrap();
        writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            alignment.query_name,
            alignment.query_length,
            alignment.query_start,
            alignment.query_end,
            alignment.strand,
            alignment.target_name,
            alignment.target_length,
            alignment.target_start,
            alignment.target_end,
            alignment.num_matching,
            alignment.alignment_block_length,
            alignment.mapq)?;
    }

    return Ok(overlaps);
}