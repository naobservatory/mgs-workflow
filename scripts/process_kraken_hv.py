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

# I/O functions

def get_virus_db(virus_path):
    """Generate data frame of human virus taxids."""
    df_viruses = pd.read_csv(virus_path, sep="\t", header=None).rename(columns={0:"taxid", 1:"name"})
    return(df_viruses)

def get_parents(nodes_path):
    """Convert a nodes DMP file into a dictionary of child:parent taxon mappings."""
    parents = {}  # child_taxid -> parent_taxid
    with open(nodes_path) as inf:
        for line in inf:
            child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
            child_taxid = int(child_taxid)
            parent_taxid = int(parent_taxid)
            parents[child_taxid] = parent_taxid
    return(parents)

# Per-line processing functions for Kraken output

def extract_assignment(name_and_taxid):
    """Extract assigned name and taxid (if any) from kraken DB field."""
    pattern = "^(.*) [(]taxid ([0-9]+)[)]$"
    name, taxid = re.findall(pattern, name_and_taxid)[0]
    return(name, int(taxid))

def screen_assignment(taxid, parents, virus_db):
    """Check whether a read has been assigned to a human-infecting virus."""
    virus_taxids = virus_db["taxid"].unique()
    while taxid not in [0,1,2]:
        if taxid in virus_taxids: return(True)
        taxid = parents[taxid]
    return(False)

def process_line(kraken_line, parents, virus_db):
    """Extract information from one line of a kraken DB and return a processed string."""
    # Strip and split line
    line = kraken_line.strip()
    if not line: return(None)
    is_classified, seq_id, name_and_taxid, length, encoded_hits = line.split("\t")
    # Process fields
    classified = True if is_classified == "C" else False if is_classified == "U" else "NA"
    assigned_name, assigned_taxid = extract_assignment(name_and_taxid)
    assigned_hv = screen_assignment(assigned_taxid, parents, virus_db)
    fields = [classified, seq_id, assigned_name, assigned_taxid, assigned_hv, length, encoded_hits]
    return(fields)

def join_line(fields):
    "Convert a list of arguments into a TSV string for output."
    return("\t".join(map(str, fields)) + "\n")

# Whole-file processing functions

def process_kraken(kraken_path, out_path, parents, virus_db):
    """Process Kraken2 output into a TSV containing additional assignment information."""
    with open_by_suffix(kraken_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Write headings
        headers = ["classified", "seq_id", "assigned_name", "assigned_taxid", "assigned_hv", "length", "encoded_hits"]
        header_line = join_line(headers)
        outf.write(header_line)
        for line in inf:
            processed_line = process_line(line, parents, virus_db)
            joined_line = join_line(processed_line)
            outf.write(joined_line)

# Main function

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process Kraken2 output into a TSV with additional information.")
    parser.add_argument("kraken_output", help="Path to Kraken2 output file.")
    parser.add_argument("virus_names", help="Path to TSV file giving taxids and names of human-infecting viruses.")
    parser.add_argument("nodes_dmp", help="Path to DMP file specifying tree structure.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    args = parser.parse_args()
    kraken_path = args.kraken_output
    virus_path = args.virus_names
    nodes_path = args.nodes_dmp
    out_path = args.output_path
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Kraken file path: {}".format(kraken_path))
    print_log("Virus DB path: {}".format(virus_path))
    print_log("Taxonomy nodes path: {}".format(nodes_path))
    print_log("Output path: {}".format(out_path))
    # Extract virus taxids
    print_log("Importing virus DB...")
    df_viruses = get_virus_db(virus_path)
    virus_taxids = list(df_viruses["taxid"])
    print_log("Virus DB imported. {} total viral taxids.".format(len(virus_taxids)))
    # Extract node mapping
    print_log("Importing taxid tree structure...")
    parents = get_parents(nodes_path)
    print_log("Tree structure imported.")
    # Process Kraken2 output
    print_log("Processing kraken2 output...")
    process_kraken(kraken_path, out_path, parents, df_viruses)
    print_log("File processed.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
