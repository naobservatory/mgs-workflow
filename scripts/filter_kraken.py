#!/usr/bin/env python

# Import modules
import re
import sys
import argparse
import pandas as pd
from collections import Counter
from multiprocessing import Pool
from functools import partial
import time
import datetime
import gzip
import bz2

def print_log(message):
    print("[", datetime.datetime.now(), "]\t", message, sep="")

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

def filter_assignment_single(kraken_line, virus_taxids):
    """Return Kraken line if assigned to a human virus, otherwise None."""
    taxid = get_assignment(kraken_line)
    if taxid in virus_taxids:
        return(kraken_line)
    else:
        return(None)

def filter_hits_single(kraken_line, virus_taxids):
    """Return Kraken line if it contains hits to a human virus, otherwise None."""
    taxids = get_hits(kraken_line)
    hv = [taxid in virus_taxids for taxid in taxids]
    if any(hv):
        return(kraken_line)
    else:
        return(None)

def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, 'r')
    else:
        return open(filename, 'r')

def filter_kraken(kraken_path, filter_fn, pool):
    """Filter Kraken2 output to reads matching human-infecting viruses."""
    lines_keep = []
    with open_by_suffix(kraken_path) as inf:
        for line in pool.map(filter_fn, inf):
            if line is not None:
                lines_keep.append(line)
    return(lines_keep)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Filter Kraken2 output to reads matching human-infecting viruses.")
    parser.add_argument("kraken_output", help="Path to Kraken2 output file.")
    parser.add_argument("virus_names", help="Path to TSV file giving taxids and names of human-infecting viruses.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    parser.add_argument("-a", "--assignments", help="Filter for reads assigned to human-infecting viruses (rather than any read containing a hit).", 
                        default = False, action = "store_true")
    parser.add_argument("-t", "--n_threads", help="Number of threads.",
                        default=1, nargs=1)
    args = parser.parse_args()
    kraken_path = args.kraken_output
    virus_path = args.virus_names
    out_path = args.output_path
    get_assignments = args.assignments
    n_threads = int(args.n_threads[0])
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Kraken file path: {}".format(kraken_path))
    print_log("Virus DB path: {}".format(virus_path))
    print_log("Output path: {}".format(out_path))
    print_log("Number of threads: {}".format(n_threads))
    screen_msg = "Screening for viral {}.".format("assignments" if get_assignments else "hits")
    print_log(screen_msg)
    # Extract virus taxids
    print_log("Importing virus DB...")
    df_viruses = get_virus_db(virus_path)
    virus_taxids = list(df_viruses["taxid"])
    print_log("Virus DB imported. {} total viral taxids.".format(len(virus_taxids)))
    # Define mapping function
    filter_fn_base = filter_assignment_single if get_assignments else filter_hits_single
    filter_fn = partial(filter_fn_base, virus_taxids = virus_taxids)
    # Filter Kraken2 output
    print_log("Processing file...")
    with Pool(n_threads) as pool:
        kraken_out = filter_kraken(kraken_path, filter_fn, pool)
    print_log("File processed. {} line(s) returned.".format(len(kraken_out)))
    # Write output
    print_log("Writing output file...")
    with open(out_path, "w") as outf:
        for line in kraken_out:
            outf.write(line)
    print_log("Output file written.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()

