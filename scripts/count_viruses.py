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

def get_counts_db(counts_path):
    """Extract data frame of Kraken2 clade counts."""
    df_counts = pd.read_csv(counts_path)
    return(df_counts)

def get_virus_counts(df_counts, df_viruses):
    """Generate labelled counts dataframe for human-infecting viruses."""
    df_counts_viral = df_viruses.merge(df_counts).sort_values("taxid")
    return(df_counts_viral)

def count_hv_reads(df_counts_viral):
    """Count all human-infecting virus reads."""
    ca = df_counts_viral["clade_assignments"].sum()
    da = df_counts_viral["direct_assignments"].sum()
    return([ca,da])

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate clade counts for human-infecting viruses.")
    parser.add_argument("counts_db", help="Path to Kraken2 clade counts CSV file.")
    parser.add_argument("virus_names", help="Path to TSV file giving taxids and names of human-infecting viruses.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    args = parser.parse_args()
    counts_path = args.counts_db
    virus_path = args.virus_names
    out_path = args.output_path
    # Get input data frames
    df_viruses = get_virus_db(virus_path)
    df_counts = get_counts_db(counts_path)
    # Generate virus counts
    df_counts_viral = get_virus_counts(df_counts, df_viruses)
    n_viral = count_hv_reads(df_counts_viral)
    print(n_viral)
    # Write output
    gzip_output = True if out_path.split(".")[-1] == "gz" else False
    if gzip_output:
        df_counts_viral.to_csv(out_path, compression = "gzip", index = False)
    else:
        df_counts_viral.to_csv(out_path, index = False)

if __name__ == "__main__":
    main()
