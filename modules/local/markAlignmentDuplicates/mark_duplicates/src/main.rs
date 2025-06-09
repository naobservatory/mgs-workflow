// ------------------------------------------------------------------------------------------------
// IMPORTS
// ------------------------------------------------------------------------------------------------

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::error::Error;
use std::cmp::Ordering;
use flate2::{Compression as GzCompression, write::GzEncoder, read::GzDecoder};
use bzip2::{Compression as BzCompression, write::BzEncoder, read::BzDecoder};
use rayon::prelude::*;
use clap::Parser;

// ------------------------------------------------------------------------------------------------
// STRUCTS AND TYPES
// ------------------------------------------------------------------------------------------------

// Minimal ReadEntry struct storing only essential data for duplicate detection
#[derive(Debug, Clone)]
struct ReadEntry {
    query_name: String,
    genome_id: String,
    aln_start: Option<i32>,
    aln_end: Option<i32>,
    avg_quality: f64,
}

// Structure to store duplicate group information without storing full read data
#[derive(Debug, Clone)]
struct DuplicateGroup {
    genome_id: String,
    exemplar_name: String,
    group_size: usize,
    pairwise_match_frac: f64,
}

// Map from query_name to (genome_id, exemplar_name) for efficient lookup during second pass
type ExemplarMap = HashMap<String, (String, String)>;

// ------------------------------------------------------------------------------------------------
// ARGUMENT PARSING
// ------------------------------------------------------------------------------------------------

/// Mark duplicate reads in alignment data
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input TSV file path
    #[arg(short, long)]
    input: String,
    /// Output database file path
    #[arg(short = 'o', long)]
    output_db: String,
    /// Output metadata file path
    #[arg(short = 'm', long)]
    output_meta: String,
    /// Position deviation tolerance (0, 1, or 2)
    #[arg(short, long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=2))]
    deviation: u8,
    /// Chunk size for parallel processing
    #[arg(short, long, default_value_t = 2000, value_parser = clap::value_parser!(u32).range(1..))]
    chunk_size: u32,
    /// Number of threads to use
    #[arg(short, long, default_value_t = 4, value_parser = clap::value_parser!(u8).range(1..))]
    num_threads: u8,
}

// ------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
// ------------------------------------------------------------------------------------------------

// Define a reader based on the file extension
fn open_reader(filename: &str) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(filename)?;
    if filename.ends_with(".gz") {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else if filename.ends_with(".bz2") {
        let decoder = BzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// Define a writer based on the file extension
fn open_writer(filename: &str) -> std::io::Result<Box<dyn Write>> {
    if filename.ends_with(".gz") {
        let file = File::create(filename)?;
        let encoder = GzEncoder::new(file, GzCompression::default());
        Ok(Box::new(BufWriter::new(encoder)))
    } else if filename.ends_with(".bz2") {
        let file = File::create(filename)?;
        let encoder = BzEncoder::new(file, BzCompression::default());
        Ok(Box::new(BufWriter::new(encoder)))
    } else {
        let file = File::create(filename)?;
        Ok(Box::new(BufWriter::new(file)))
    }
}

// Implement a custom match function for comparing ReadEntries
// (Not a valid equality relation as not transitive)
fn match_reads(a: &ReadEntry, b: &ReadEntry) -> bool {
    a.genome_id == b.genome_id &&
    unsafe { compare_positions(a.aln_start, b.aln_start, DEVIATION) } &&
    unsafe { compare_positions(a.aln_end, b.aln_end, DEVIATION) }
}

// Compare the positions with a deviation
fn compare_positions(a: Option<i32>, b: Option<i32>, deviation: u8) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => (x - y).abs() <= deviation as i32,
        _ => false,
    }
}

// Implement ordered comparison for ReadEntry
fn compare_reads(a: &ReadEntry, b: &ReadEntry) -> Ordering {
    // Compare by average quality score
    let quality_cmp = a.avg_quality.partial_cmp(&b.avg_quality).unwrap_or(Ordering::Equal);
    // If equal, compare by query name (in reverse order)
    if quality_cmp == Ordering::Equal {
        b.query_name.cmp(&a.query_name)
    } else {
        quality_cmp
    }
}

// Parse the integer value or return None if the value is "NA"
fn parse_int_or_na(s: &str) -> Option<i32> {
    if s == "NA" {
        None
    } else {
        s.parse().ok()
    }
}

// Convert the ASCII quality score to a quality score (optimized for speed)
fn ascii_to_quality_score(ascii_score: &str) -> f64 {
    let bytes = ascii_score.as_bytes();
    let sum: u32 = bytes.iter().map(|&b| (b - 33) as u32).sum();
    sum as f64 / bytes.len() as f64
}

// Calculate the average quality score of the forward and reverse reads
fn average_quality_score(quality_fwd: &str, quality_rev: &str) -> f64 {
    let fwd_score = ascii_to_quality_score(quality_fwd);
    let rev_score = ascii_to_quality_score(quality_rev);
    (fwd_score + rev_score) / 2.0
}

// ------------------------------------------------------------------------------------------------
// EXTRACTION FUNCTIONS
// ------------------------------------------------------------------------------------------------

fn process_header_line(line: &str) -> Result<(Vec<&str>, HashMap<&str, usize>, usize), Box<dyn Error>> {
    // Split the line by tabs and collect the headers
    let headers: Vec<&str> = line.split('\t').collect();
    let header_count: usize = headers.len();
    // Build a map from header fields to indices
    let header_indices: HashMap<_, _> = headers.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    // Define required header fields
    let required_headers = vec![
        "seq_id", "bowtie2_genome_id_all", "bowtie2_ref_start_fwd", "bowtie2_ref_start_rev",
        "query_qual_fwd", "query_qual_rev"
    ];
    // Build a lookup for required headers
    let mut indices = HashMap::new();
    for header in required_headers {
        let idx = header_indices.get(header)
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Missing required header: {}", header)))?;
        indices.insert(header, *idx);
    }
    // Return output
    Ok((headers, indices, header_count))
}

// Efficient function that creates ReadEntry with minimal memory allocation
fn make_read_entry(fields: &[String], indices: &HashMap<&str, usize>) -> ReadEntry {
    // Extract required fields using references to avoid cloning unnecessarily
    let query_name = fields[indices["seq_id"]].clone();
    let genome_id = &fields[indices["bowtie2_genome_id_all"]];
    let ref_start_fwd = parse_int_or_na(&fields[indices["bowtie2_ref_start_fwd"]]);
    let ref_start_rev = parse_int_or_na(&fields[indices["bowtie2_ref_start_rev"]]);
    let quality_fwd = &fields[indices["query_qual_fwd"]];
    let quality_rev = &fields[indices["query_qual_rev"]];
    // Handle split assignments
    let genome_id_sorted: String;
    let aln_start: Option<i32>;
    let aln_end: Option<i32>;
    if genome_id.contains('/') {
        // Split genome_id by "/", sort the parts, and join them
        let parts: Vec<&str> = genome_id.split('/').collect();
        let mut sorted_parts = parts.clone();
        sorted_parts.sort();
        genome_id_sorted = sorted_parts.join("/");
        // Get the index of the first genome ID in the sorted list
        let genome_id_index = sorted_parts.iter().position(|&s| s == parts[0]).unwrap();
        // Arrange start coordinates to correspond to sorted genome IDs
        (aln_start, aln_end) = if genome_id_index == 0 {
            (ref_start_fwd, ref_start_rev)
        } else {
            (ref_start_rev, ref_start_fwd)
        };
    } else {
        // If only one genome ID, use it directly
        genome_id_sorted = genome_id.to_string();
        // Normalize coordinates: if values are present, use the minimum and maximum
        (aln_start, aln_end) = match (ref_start_fwd, ref_start_rev) {
            (Some(fwd), Some(rev)) => (Some(fwd.min(rev)), Some(fwd.max(rev))),
            _ => (ref_start_fwd, ref_start_rev),
        };
    };
    let avg_quality = average_quality_score(quality_fwd, quality_rev);
    // Return the ReadEntry with minimal memory footprint
    ReadEntry { 
        query_name, 
        genome_id: genome_id_sorted, 
        aln_start, 
        aln_end, 
        avg_quality 
    }
}

// Process a chunk of lines in parallel to create ReadEntry objects
fn process_chunk_parallel(
    lines: &[String], 
    indices: &HashMap<&str, usize>,
    header_count: usize
) -> Result<Vec<ReadEntry>, Box<dyn Error>> {
    // Parse lines in parallel using rayon
    let read_entries: Result<Vec<ReadEntry>, String> = lines
        .par_iter()  // Parallel iterator from rayon
        .map(|line| {
            // Split line into fields
            let fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            // Validate field count
            if fields.len() != header_count {
                return Err(format!("Invalid field count: {} (expected {})", fields.len(), header_count));
            }
            // Create ReadEntry from fields
            Ok(make_read_entry(&fields, indices))
        })
        .collect();
    // Convert String errors to Box<dyn Error>
    read_entries.map_err(|e| -> Box<dyn Error> { 
        std::io::Error::new(std::io::ErrorKind::InvalidData, e).into() 
    })
}

fn add_read_to_groups(
    mut groups: HashMap<String, Vec<Vec<ReadEntry>>>,
    read_entry: ReadEntry,
) -> HashMap<String, Vec<Vec<ReadEntry>>> {
    // Check if the genome_id is already present
    if !groups.contains_key(&read_entry.genome_id) {
        // If not, create a new ID entry and add the read entry to it as a new group
        groups.insert(read_entry.genome_id.clone(), Vec::new());
        groups.get_mut(&read_entry.genome_id).unwrap().push(vec![read_entry]);
    } else {
        // Otherwise, compare the read entry with each group in the ID entry
        let mut groups_match = Vec::new();
        for i in 0.. groups[&read_entry.genome_id].len() {
            for j in 0..groups[&read_entry.genome_id][i].len() {
                if match_reads(&read_entry, &groups[&read_entry.genome_id][i][j]) {
                    groups_match.push(i);
                    break;
                }
            }
        }
        // If no matches, create a new group in the ID entry and add the read entry to it
        if groups_match.is_empty() {
            groups.get_mut(&read_entry.genome_id).unwrap().push(vec![read_entry]);
        // If exactly one match, add the read entry to that group
        } else if groups_match.len() == 1 {
            groups.get_mut(&read_entry.genome_id).unwrap()[groups_match[0]].push(read_entry);
        // If multiple matches, combine them into a new group, add the read entry, then
        // delete the old groups
        } else {
            // Combine the groups into a new group
            let mut combined_group = Vec::new();
            for i in groups_match.iter() {
                combined_group.extend(groups[&read_entry.genome_id][*i].clone());
            }
            // Add the read entry to the combined group
            combined_group.push(read_entry.clone());
            // Add the combined group to the ID entry
            groups.get_mut(&read_entry.genome_id).unwrap().push(combined_group);
            // Remove the old groups in reverse order (to avoid index shifting)
            for i in groups_match.iter().rev() {
                groups.get_mut(&read_entry.genome_id).unwrap().remove(*i);
            }
        }
    }
    // Return the updated groups
    groups
}

fn extract_read_groups(input_path: &str,
    chunk_size: u32
) -> Result<(String, HashMap<String, Vec<Vec<ReadEntry>>>, usize), Box<dyn Error>> {
    // Open the input file
    let reader = open_reader(input_path)?;
    // Create a HashMap to store duplicate information
    let mut groups: HashMap<String, Vec<Vec<ReadEntry>>> = HashMap::new();
    // Process the header line and derive the required fields
    let mut lines = reader.lines();
    let header_line = lines.next().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty input file"))??;
    let (headers, indices, header_count) = process_header_line(&header_line)?;
    // Create the output header line
    let mut headers_out = headers.clone();
    headers_out.push("bowtie2_dup_exemplar");
    let header_out = headers_out.join("\t");
    // Get the seq_id column index for later use
    let seq_id_index = indices["seq_id"];
    // Read and process the input file in chunks
    let mut line_buffer = Vec::new();
    for line in lines {
        let line = line?;
        line_buffer.push(line);
        // Process chunk when buffer is full
        if line_buffer.len() >= chunk_size as usize {
            // Process this chunk in parallel
            let read_entries = process_chunk_parallel(&line_buffer, &indices, header_count)?;
            // Add to groups sequentially (to avoid race conditions)
            for read_entry in read_entries {
                groups = add_read_to_groups(groups, read_entry);
            }
            // Clear the buffer
            line_buffer.clear();
        }
    }
    // Process remaining lines in the buffer
    if !line_buffer.is_empty() {
        let read_entries = process_chunk_parallel(&line_buffer, &indices, header_count)?;
        for read_entry in read_entries {
            groups = add_read_to_groups(groups, read_entry);
        }
    }
    Ok((header_out, groups, seq_id_index))
}

// ------------------------------------------------------------------------------------------------
// PROCESSING FUNCTIONS
// ------------------------------------------------------------------------------------------------

// Process duplicate groups to create exemplar mapping and metadata (focused on group processing)
fn process_read_groups(
    groups: HashMap<String, Vec<Vec<ReadEntry>>>
) -> Result<(ExemplarMap, Vec<DuplicateGroup>), Box<dyn Error>> {
    let mut exemplar_map = ExemplarMap::new();
    let mut duplicate_groups = Vec::new();
    for (genome_id, id_group) in groups {
        for dup_group in id_group {
            // Find the exemplar using compare_reads
            let exemplar = dup_group.iter().max_by(|a, b| compare_reads(a, b)).unwrap();
            let exemplar_name = exemplar.query_name.clone();
            // Calculate size of duplicate group
            let dup_count = dup_group.len();
            // Calculate fraction of pairwise matches (as a QC metric for the group as a whole)
            let pairwise_match_frac: f64;
            if dup_count == 1 {
                pairwise_match_frac = 1.0;
            } else {
                let mut pairwise_match_count: f64 = 0.0;
                let dup_count_float: f64 = dup_count as f64;
                let n_pairs: f64 = dup_count_float * (dup_count_float - 1.0) / 2.0;
                for i in 0..dup_count {
                    for j in (i + 1)..dup_count {
                        let read_i = &dup_group[i];
                        let read_j = &dup_group[j];
                        if match_reads(read_i, read_j) {
                            pairwise_match_count += 1.0;
                        }
                    }
                }
                pairwise_match_frac = pairwise_match_count / n_pairs;
            }
            // Store duplicate group information
            duplicate_groups.push(DuplicateGroup {
                genome_id: genome_id.clone(),
                exemplar_name: exemplar_name.clone(),
                group_size: dup_count,
                pairwise_match_frac,
            });
            // Map each read in the group to its exemplar
            for read_entry in dup_group {
                exemplar_map.insert(read_entry.query_name, (genome_id.clone(), exemplar_name.clone()));
            }
        }
    }
    Ok((exemplar_map, duplicate_groups))
}

// ------------------------------------------------------------------------------------------------
// WRITING FUNCTIONS
// ------------------------------------------------------------------------------------------------

// Write duplicate group metadata file (no file streaming required)
fn write_metadata_file(
    duplicate_groups: &Vec<DuplicateGroup>,
    output_path_meta: &str,
) -> Result<(), Box<dyn Error>> {
    // Open the metadata output file
    let mut writer_meta = open_writer(output_path_meta)?;
    // Write header
    let header_meta = "bowtie2_genome_id_all\tbowtie2_dup_exemplar\tbowtie2_dup_count\tbowtie2_dup_pairwise_match_frac";
    writeln!(writer_meta, "{}", header_meta)?;
    // Write duplicate group metadata (once per group)
    for dup_group in duplicate_groups {
        writeln!(writer_meta, "{}\t{}\t{}\t{}", 
                dup_group.genome_id, dup_group.exemplar_name, dup_group.group_size, dup_group.pairwise_match_frac)?;
    }
    Ok(())
}

// Stream through file and add exemplar information
fn write_database_file(
    input_path: &str,
    header_out: &str,
    exemplar_map: &ExemplarMap,
    seq_id_index: usize,
    output_path_db: &str,
) -> Result<(), Box<dyn Error>> {
    // Open input file for second pass
    let reader = open_reader(input_path)?;
    // Open the database output file
    let mut writer_db = open_writer(output_path_db)?;
    // Write header
    writeln!(writer_db, "{}", header_out)?;
    // Process input file line by line for output generation
    let mut lines = reader.lines();
    let _header_line = lines.next(); // Skip header
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let query_name = fields[seq_id_index];
        // Look up exemplar for this read
        if let Some((_genome_id, exemplar_name)) = exemplar_map.get(query_name) {
            writeln!(writer_db, "{}\t{}", line, exemplar_name)?;
        } else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Could not find exemplar for read: {}", query_name)
            ).into());
        }
    }
    Ok(())
}

// ------------------------------------------------------------------------------------------------
// TOP-LEVEL FUNCTIONS
// ------------------------------------------------------------------------------------------------

// Define the deviation value
static mut DEVIATION: u8 = 0;

// Two-pass processing for improved memory efficiency
fn process_tsv(input_path: &str,
    output_path_db: &str,
    output_path_meta: &str,
    chunk_size: u32) -> Result<(), Box<dyn Error>> {
    // Extract read groups from the input file
    let (header_out, groups, seq_id_index) = extract_read_groups(input_path, chunk_size)?;
    // Process duplicate groups to create exemplar mapping and metadata
    let (exemplar_map, duplicate_groups) = process_read_groups(groups)?;
    // Write metadata file
    write_metadata_file(&duplicate_groups, output_path_meta)?;
    // Write database file
    write_database_file(input_path, &header_out, &exemplar_map, seq_id_index, output_path_db)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command line arguments
    let args = Args::parse();
    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads as usize)
        .build_global()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, 
            format!("Failed to configure thread pool: {}", e)))?;
    // Set the deviation value
    unsafe {
        DEVIATION = args.deviation;
    }
    // Run the main processing function
    return process_tsv(&args.input, &args.output_db, &args.output_meta, args.chunk_size);
}
