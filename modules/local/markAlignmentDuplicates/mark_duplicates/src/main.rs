use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::cmp::Ordering;
use flate2::{Compression as GzCompression, write::GzEncoder, read::GzDecoder};
use bzip2::{Compression as BzCompression, write::BzEncoder, read::BzDecoder};

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

// Define a ReadEntry struct to store the read information
#[derive(Debug, Clone)]
struct ReadEntry {
    query_name: String,
    genome_id: String,
    aln_start: Option<i32>,
    aln_end: Option<i32>,
    avg_quality: f64,
    fields: Vec<String>,
}

// Implement a custom match function for comparing ReadEntries
// (Not a valid equality relation as not transitive)
fn match_reads(a: &ReadEntry, b: &ReadEntry) -> bool {
    a.genome_id == b.genome_id &&
    unsafe { compare_positions(a.aln_start, b.aln_start, DEVIATION) } &&
    unsafe { compare_positions(a.aln_end, b.aln_end, DEVIATION) }
}

// Compare the positions with a deviation
// If both positions are None, they are considered equal
fn compare_positions(a: Option<i32>, b: Option<i32>, deviation: i32) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => (x - y).abs() <= deviation,
        (None, None) => true,
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

// Convert the ASCII quality score to a quality score
fn ascii_to_quality_score(ascii_score: &str) -> f64 {
    if ascii_score == "NA" {
        return 0.0;
    }
    ascii_score.chars()
        .map(|c| c as u8 as f64 - 33.0)
        .sum::<f64>() / ascii_score.len() as f64
}

// Calculate the average quality score of the forward and reverse reads
fn average_quality_score(quality_fwd: &str, quality_rev: &str) -> f64 {
    let fwd_score = ascii_to_quality_score(quality_fwd);
    let rev_score = ascii_to_quality_score(quality_rev);
    (fwd_score + rev_score) / 2.0
}

fn process_header_line(line: &str) -> Result<(Vec<&str>, HashMap<&str, usize>, usize), Box<dyn Error>> {
    // Split the line by tabs and collect the headers
    let headers: Vec<&str> = line.split('\t').collect();
    let header_count: usize = headers.len();
    // Build a map from header fields to indices
    let header_indices: HashMap<_, _> = headers.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    // Define required header fields
    let required_headers = vec![
        "seq_id", "aligner_genome_id_all", "aligner_ref_start", "aligner_ref_start_rev",
        "query_qual", "query_qual_rev"
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

fn make_read_entry(fields: Vec<String>, indices: &HashMap<&str, usize>) -> ReadEntry {
    // Extract required fields
    let query_name = fields[indices["seq_id"]].to_string();
    let genome_id = fields[indices["aligner_genome_id_all"]].to_string();
    let ref_start_fwd = parse_int_or_na(&fields[indices["aligner_ref_start"]]);
    let ref_start_rev = parse_int_or_na(&fields[indices["aligner_ref_start_rev"]]);
    let quality_fwd = fields[indices["query_qual"]].to_string();
    let quality_rev = fields[indices["query_qual_rev"]].to_string();
    
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
        // Note: this doesn't need to handle the case where one value is None because then you could never get multiple genome_ids
        (aln_start, aln_end) = if genome_id_index == 0 {
            (ref_start_fwd, ref_start_rev)
        } else {
            (ref_start_rev, ref_start_fwd)
        };
    } else {
        // If only one genome ID, use it directly
        genome_id_sorted = genome_id;
        // Normalize coordinates: handle cases where one value is None
        (aln_start, aln_end) = match (ref_start_fwd, ref_start_rev) {
            (Some(fwd), Some(rev)) => (Some(fwd.min(rev)), Some(fwd.max(rev))),
            (Some(fwd), None) => (Some(fwd), None),
            (None, Some(rev)) => (Some(rev), None),
            (None, None) => (None, None),
        };
    };
    let avg_quality = average_quality_score(&quality_fwd, &quality_rev);
    // Return the ReadEntry
    ReadEntry { query_name, genome_id: genome_id_sorted, aln_start, aln_end, avg_quality, fields }
}

fn extract_read_groups(input_path: &str) -> Result<(String, HashMap<String, Vec<Vec<ReadEntry>>>), Box<dyn Error>> {
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
    // Read and process the input file
    for (index, line) in lines.enumerate() {
        // Extract fields and verify count
        let line = line?;
        let fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if fields.len() != header_count {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Incorrect field count in line {}: {} (observed) vs {} (expected)", 
                    index + 1, fields.len(), header_count),
            ).into());
        }
        // Create a ReadEntry from the fields
        let read_entry = make_read_entry(fields.clone(), &indices);
        // Add the read to the groups HashMap
        groups = add_read_to_groups(groups, read_entry);
    }
    Ok((header_out, groups))
}

fn process_read_groups(header_out: &str, groups: HashMap<String, Vec<Vec<ReadEntry>>>,
    output_path_db: &str, output_path_meta: &str) -> Result<(), Box<dyn Error>> {
    // Open the output files
    let mut writer_db = open_writer(output_path_db)?;
    let mut writer_meta = open_writer(output_path_meta)?;
    // Write the header line to the output files
    let header_meta = "aligner_genome_id_all\tbowtie2_dup_exemplar\tbowtie2_dup_count\tbowtie2_dup_pairwise_match_frac";
    writeln!(writer_db, "{}", header_out)?;
    writeln!(writer_meta, "{}", header_meta)?;
    // Process each group of reads
    for (id, id_group) in groups {
        for dup_group in id_group {
            // Find the exemplar using compare_reads
            let exemplar = dup_group.iter().max_by(|a, b| compare_reads(a, b)).unwrap().query_name.clone();
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
                        // Match using ReadEntry PartialEq implementation
                        if match_reads(read_i, read_j) {
                            pairwise_match_count += 1.0;
                        }
                    }
                }
                pairwise_match_frac = pairwise_match_count / n_pairs;
            }
            // Write duplicate group entry to metadata file
            writeln!(writer_meta, "{}\t{}\t{}\t{}", id, exemplar, dup_count, pairwise_match_frac)?;
            // Write individual entries to main output file
            for read_entry in dup_group {
                writeln!(writer_db, "{}\t{}", read_entry.fields.join("\t"), exemplar)?;
            }
        }
    }
    Ok(())
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

fn process_tsv(input_path: &str, output_path_db: &str, output_path_meta: &str) -> Result<(), Box<dyn Error>> {
    // Extract read groups from the input file
    let (header_out, groups) = extract_read_groups(input_path)?;
    // Process read groups and write output
    process_read_groups(&header_out, groups, output_path_db, output_path_meta)?;
    Ok(())
}

// Define the deviation value
static mut DEVIATION: i32 = 0;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    // Check if the correct number of arguments are provided
    if args.len() != 5 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!(
            "Usage: {} <input_file> <output_file_db> <output_file_meta> <deviation>", args[0]
        )).into());
    }
    // Parse the input and output file paths
    let input_path = &args[1];
    let output_path_db = &args[2];
    let output_path_meta = &args[3];
    // Parse the deviation value
    let deviation: i32 = args[4].parse().map_err(|_| {
        std::io::Error::new(std::io::ErrorKind::InvalidInput, "Deviation must be an integer")
    })?;
    // Check if the deviation value is valid
    if deviation != 0 && deviation != 1 && deviation != 2 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!(
            "Deviation must be at most 2"
        )).into());
    }
    // Set the deviation value
    unsafe {
        DEVIATION = deviation;
    }
    // Run the main processing function
    return process_tsv(input_path, output_path_db, output_path_meta);
}
