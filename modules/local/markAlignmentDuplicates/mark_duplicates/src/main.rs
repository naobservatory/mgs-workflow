use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::env;
use std::hash::{Hash, Hasher};
use std::error::Error;
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

// Defune a writer based on the file extension
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

// Define the PositionKey struct to store the genome_id, forward and reverse reference start positions
#[derive(Debug, Clone, Eq)]
struct PositionKey {
    genome_id: String,
    aln_start: Option<i32>,
    aln_end: Option<i32>,
}

// Implement the PartialEq trait for PositionKey
impl PartialEq for PositionKey {
    fn eq(&self, other: &Self) -> bool {
        self.genome_id == other.genome_id &&
        unsafe { compare_positions(self.aln_start, other.aln_start, DEVIATION) } &&
        unsafe { compare_positions(self.aln_end, other.aln_end, DEVIATION) }
    }
}

// Compare the positions with a deviation
fn compare_positions(a: Option<i32>, b: Option<i32>, deviation: i32) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => (x - y).abs() <= deviation,
        _ => false,
    }
}

// Implement the Hash trait for PositionKey
impl Hash for PositionKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.genome_id.hash(state);
        unsafe {
            hash_position(self.aln_start, DEVIATION, state);
            hash_position(self.aln_end, DEVIATION, state);
        }
    }
}

// Hash the position with a deviation
fn hash_position<H: Hasher>(pos: Option<i32>, deviation: i32, state: &mut H) {
    if let Some(p) = pos {
        ((p + deviation) / (2 * deviation + 1)).hash(state);
    } else {
        None::<i32>.hash(state);
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

fn process_tsv(input_path: &str, output_path: &str) -> Result<(), Box<dyn Error>> {
    // Open the input file
    let reader = open_reader(input_path)?;
    let mut writer = open_writer(output_path)?;
    // Create a HashMap to store duplicate information
    let mut duplicates: HashMap<PositionKey, Vec<(String, f64, Vec<String>)>> = HashMap::new();
    // Skip the header line
    let mut lines = reader.lines();
    let header_line = lines.next().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty input file"))??;
    let headers: Vec<&str> = header_line.split('\t').collect();
    let header_count = headers.len(); // Number of header fields
    // Build a map from header fields to indices
    let header_indices: std::collections::HashMap<_, _> = headers.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    // Define required header fields
    let required_headers = vec![
        "seq_id", "bowtie2_genome_id_all", "bowtie2_ref_start_fwd", "bowtie2_ref_start_rev",
        "query_qual_fwd", "query_qual_rev"
    ];
    // Build a lookup for required headers
    let mut indices = std::collections::HashMap::new();
    for header in required_headers {
        let idx = header_indices.get(header)
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Missing required header: {}", header)))?;
        indices.insert(header, *idx);
    }
    // Write header with additional columns
    writeln!(writer, "{}\t{}\t{}", headers.join("\t"), "bowtie2_dup_exemplar", "bowtie2_dup_count")?;
    // Read and process the input file
    for (index, line) in lines.enumerate() {
        let line = line?;
        let fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if fields.len() != header_count {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Incorrect field count in line {}: {} (observed) vs {} (expected)", index + 1, fields.len(), header_count),
            ).into());
        }
        let query_name = fields[indices["seq_id"]].to_string();
        let genome_id = fields[indices["bowtie2_genome_id_all"]].to_string();
        let ref_start_fwd = parse_int_or_na(&fields[indices["bowtie2_ref_start_fwd"]]);
        let ref_start_rev = parse_int_or_na(&fields[indices["bowtie2_ref_start_rev"]]);
        let quality_fwd = fields[indices["query_qual_fwd"]].to_string();
        let quality_rev = fields[indices["query_qual_rev"]].to_string();
        // Normalize the coordinates: if both values are present, use the minimum and maximum
        let (aln_start, aln_end) = match (ref_start_fwd, ref_start_rev) {
            (Some(fwd), Some(rev)) => (Some(fwd.min(rev)), Some(fwd.max(rev))),
            _ => (ref_start_fwd, ref_start_rev),
        };
        // Calculate the average quality score of the forward and reverse reads
        let avg_quality = average_quality_score(&quality_fwd, &quality_rev);
        // Create a PositionKey for the current read
        let key = PositionKey {
            genome_id: genome_id.clone(),
            aln_start: aln_start,
            aln_end: aln_end,
        };
        // Add the read to the duplicates HashMap
        duplicates.entry(key).or_default().push((query_name, avg_quality, fields));
    }
    // Process duplicates and write output
    for (_key, group) in duplicates {
        // Find the exemplar as the read with the highest average quality score
        // If reads tie, take the lexicographically first one
        let exemplar = group.iter()
            .max_by(|a, b| {
                let quality_cmp = a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal);
                if quality_cmp == std::cmp::Ordering::Equal { b.0.cmp(&a.0) } else { quality_cmp }
            })
            .unwrap()
            .0.clone();
        let dup_count = group.len();
        for (_query_name, _avg_quality, fields) in group {
            writeln!(
                writer,
                "{}\t{}\t{}",
                fields.join("\t"),
                exemplar,
                dup_count
            )?;
        }
    }
    Ok(())
}

// Define the deviation value
static mut DEVIATION: i32 = 0;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    // Check if the correct number of arguments are provided
    if args.len() != 4 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!(
            "Usage: {} <input_file> <output_file> <deviation>", args[0]
        )).into());
    }
    // Parse the input and output file paths
    let input_path = &args[1];
    let output_path = &args[2];
    // Parse the deviation value
    let deviation: i32 = args[3].parse().map_err(|_| {
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
    return process_tsv(input_path, output_path)
}
