#!/usr/bin/env python

# Import modules
import re
import sys
import argparse
import pandas as pd
from collections import Counter

def get_virus_db(virus_path):
    """Generate data frame of human virus taxids."""
    df_viruses = pd.read_csv(virus_path, sep="\t", header=None).rename(columns={0:"taxid", 1:"name"})
    return(df_viruses)

def get_assignment(kraken_line):
    """Extract the assigned taxid from one line of a kraken DB."""
    line = kraken_line.strip()
    if not line: return(None)
    _, _, name_and_taxid, _, encoded_hits = line.split("\t")
    taxid, = re.findall("^.*[(]taxid ([0-9]+)[)]$", name_and_taxid)
    taxid = int(taxid)
    return(taxid)

def get_hits(kraken_line):
    """Extract taxid hits from one line of a kraken DB."""
    line = kraken_line.strip()
    if not line: return(None)
    taxids = []
    try:
        taxid_matches = line.strip().split("\t")[4]
        for taxid_match in taxid_matches.split(" "):
            taxid, n_kmers = taxid_match.split(":")
            if taxid == "A": continue # Ambiguous nucleotide
            if taxid == "|": continue # Paired end transition
            taxid = int(taxid)
            taxids.append(taxid)
    except Exception:
        print(line)
        raise
    return(taxids)

def filter_kraken_hits(kraken_path, virus_db):
    """Filter Kraken2 output to reads containing hits to human-infecting viruses."""
    human_viruses = list(virus_db["taxid"])
    lines_keep = []
    with open(kraken_path) as inf:
        for line in inf:
            keep = False
            taxids = get_hits(line)
            for taxid in taxids:
                if taxid in human_viruses:
                    keep = True
            if keep:
                lines_keep.append(line)
    return(lines_keep)

def filter_kraken_assignments(kraken_path, virus_db):
    """Filter Kraken2 output to reads assigned to human-infecting viruses."""
    human_viruses = list(virus_db["taxid"])
    lines_keep = []
    with open(kraken_path) as inf:
        for line in inf:
            taxid = get_assignment(line)
            if taxid in human_viruses:
                lines_keep.append(line)
    return(lines_keep)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Filter Kraken2 output to reads matching human-infecting viruses.")
    parser.add_argument("kraken_output", help="Path to Kraken2 output file.")
    parser.add_argument("virus_names", help="Path to TSV file giving taxids and names of human-infecting viruses.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    parser.add_argument("--assignments", help="Filter for reads assigned to human-infecting viruses (rather than any read containing a hit).", 
                        default = False, action = "store_true")
    args = parser.parse_args()
    kraken_path = args.kraken_output
    virus_path = args.virus_names
    out_path = args.output_path
    get_assignments = args.assignments
    # Extract virus taxids
    df_viruses = get_virus_db(virus_path)
    # Filter Kraken2 output
    filter_fn = filter_kraken_assignments if get_assignments else filter_kraken_hits
    kraken_out = filter_fn(kraken_path, df_viruses)
    # Write output
    with open(out_path, "w") as outf:
        for line in kraken_out:
            outf.write(line)

if __name__ == "__main__":
    main()

