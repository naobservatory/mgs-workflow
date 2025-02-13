use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::env;
use std::hash::{Hash, Hasher};

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
        //(None, None) => true,
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

fn process_tsv(input_path: &str, output_path: &str) -> std::io::Result<()> {
    // Open the input file
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);
    let mut writer = BufWriter::new(File::create(output_path)?);
    // Create a HashMap to store duplicate information
    let mut duplicates: HashMap<PositionKey, Vec<(String, Option<i32>, f64)>> = HashMap::new();
    // Skip the header line
    let mut lines = reader.lines();
    let _header = lines.next().ok_or(std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty file"))??;
    // Write header (corrected)
    writeln!(writer, "query_name\tfragment_length\texemplar\tdup_count")?;
     // Read and process the input file
     for (index, line) in lines.enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        // Check if the line has the correct number of fields
        if fields.len() < 23 {
            eprintln!("Line {} has fewer than 23 fields: {:?}", index + 1, fields);
            continue;
        }
        let query_name = fields[0].to_string();
        let genome_id = fields[1].to_string();
        let fragment_length = parse_int_or_na(fields[3]);
        let aln_start = parse_int_or_na(fields[10]);
        let aln_end = parse_int_or_na(fields[11]);
        let quality_fwd = fields[20].to_string();
        let quality_rev = fields[21].to_string();
        // Calculate the average quality score of the forward and reverse reads
        let avg_quality = average_quality_score(&quality_fwd, &quality_rev);
        // Create a PositionKey for the current read
        let key = PositionKey {
            genome_id: genome_id.clone(),
            aln_start,
            aln_end,
        };
        // Add the read to the duplicates HashMap
        duplicates.entry(key).or_default().push((query_name, fragment_length, avg_quality));
    }
    // Process duplicates and write output
    for (_key, group) in duplicates {
        // Find the exemplar (read with the highest average quality score). 4 refers to the index of the avg_quality in the tuple.
        let exemplar = group.iter()
            .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap()
            .0.clone();
        let dup_count = group.len();
        for (query_name, fragment_length, _avg_quality) in group {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                query_name,
                fragment_length.unwrap_or(-1), // Use -1 for NA values
                exemplar,
                dup_count
            )?;
        }
    }
    Ok(())
}

// Define the deviation value
static mut DEVIATION: i32 = 0;

fn main() -> std::io::Result<()> {
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
            "Error: deviation must be 0, 1, or 2"
        )).into());
    }
    // Set the deviation value
    unsafe {
        DEVIATION = deviation;
    }
    // Run the main processing function
    match process_tsv(input_path, output_path) {
        Ok(_) => println!("Processing complete."),
        Err(e) => eprintln!("Error processing file: {:?}", e),
    }
    Ok(())
}
