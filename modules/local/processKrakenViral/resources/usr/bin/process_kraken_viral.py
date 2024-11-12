#!/usr/bin/env python

# Import modules
import re
import argparse
import pandas as pd
import time
import datetime
import gzip
import bz2

# Utility functions

def print_log(message):
    print("[", datetime.datetime.now(), "]\t", message, sep="")

def open_by_suffix(filename, mode="r"):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

# Per-line processing functions for Kraken output

def extract_assignment(name_and_taxid):
    """Extract assigned name and taxid (if any) from kraken DB field."""
    pattern = "^(.*) [(]taxid ([0-9]+)[)]$"
    name, taxid = re.findall(pattern, name_and_taxid)[0]
    return(name, taxid)

def process_line(kraken_line, virus_status_dict):
    """Extract information from one line of a kraken DB and return a processed string."""
    # Strip and split line
    line = kraken_line.strip()
    if not line: return(None)
    is_classified, seq_id, name_and_taxid, length, encoded_hits = line.split("\t")
    # Process fields
    classified = True if is_classified == "C" else False if is_classified == "U" else "NA"
    assigned_name, assigned_taxid = extract_assignment(name_and_taxid)
    try:
        assigned_host_virus = virus_status_dict[assigned_taxid]
    except KeyError:
        assigned_host_virus = "0"
    fields = [classified, seq_id, assigned_name, assigned_taxid, assigned_host_virus, length, encoded_hits]
    return(fields)

def join_line(fields):
    "Convert a list of arguments into a TSV string for output."
    return("\t".join(map(str, fields)) + "\n")

# Whole-file processing functions

def process_kraken(kraken_path, out_path, virus_status_dict):
    """Process Kraken2 output into a TSV containing additional assignment information."""
    with open_by_suffix(kraken_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Write headings
        headers = ["kraken_classified", "seq_id", "kraken_assigned_name", "kraken_assigned_taxid", "kraken_assigned_host_virus", "kraken_length", "kraken_encoded_hits"]
        header_line = join_line(headers)
        outf.write(header_line)
        for line in inf:
            processed_line = process_line(line, virus_status_dict)
            joined_line = join_line(processed_line)
            outf.write(joined_line)

# Main function

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process Kraken2 output into a TSV with additional information.")
    parser.add_argument("kraken_output", help="Path to Kraken2 output file.")
    parser.add_argument("virus_db", help="Path to TSV file giving viral taxonomic information.")
    parser.add_argument("host_taxon", help="Virus host taxon to screen against (e.g. human, vertebrate).")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    args = parser.parse_args()
    kraken_path = args.kraken_output
    virus_path = args.virus_db
    host_taxon = args.host_taxon
    out_path = args.output_path
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Kraken file path: {}".format(kraken_path))
    print_log("Virus DB path: {}".format(virus_path))
    print_log("Host taxon: {}".format(host_taxon))
    print_log("Output path: {}".format(out_path))
    host_column = "infection_status_" + host_taxon
    print_log("Host taxon search column: {}".format(host_column))
    # Extract virus taxids
    print_log("Importing viral DB file...")
    virus_db = pd.read_csv(virus_path, sep="\t", dtype=str)
    print_log("Virus DB imported. {} total viral taxids.".format(len(virus_db)))
    virus_status_dict = {taxid: status for taxid, status in
                         zip(virus_db["taxid"], virus_db[host_column])}
    # Process Kraken2 output
    print_log("Processing kraken2 output...")
    process_kraken(kraken_path, out_path, virus_status_dict)
    print_log("File processed.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
