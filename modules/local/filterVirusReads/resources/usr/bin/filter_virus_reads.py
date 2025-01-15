#!/usr/bin/env python

# Import modules
import argparse
import time
import datetime
import gzip
import bz2
import os

def print_log(message):
    print("[", datetime.datetime.now(), "]  ", message, sep="")

def open_by_suffix(filename, mode="r", debug=False):
    if debug:
        print_log(f"\tOpening file object: {filename}")
        print_log(f"\tOpening mode: {mode}")
        print_log(f"\tGZIP mode: {filename.endswith('.gz')}")
        print_log(f"\tBZ2 mode: {filename.endswith('.bz2')}")
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

def filter_virus_reads(input_path, score_threshold, out_path):
    """Filter putative virus reads based on Kraken and Bowtie2 results."""
    with open_by_suffix(input_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Read and handle header line
        headers = inf.readline().strip().split("\t")
        idx = {h: i for i, h in enumerate(headers)}
        # Check for necessary columns in header
        if "kraken_classified" not in headers:
            raise ValueError("Missing column in input TSV: 'kraken_classified'")
        if "kraken_assigned_host_virus" not in headers:
            raise ValueError("Missing column in input TSV: 'kraken_assigned_host_virus'")
        if "bowtie2_length_normalized_score_max" not in headers:
            raise ValueError("Missing column in input TSV: 'bowtie2_length_normalized_score_max'")
        outf.write("\t".join(headers) + "\n")
        # Read in lines and filter based on Kraken assignment and score
        for line in inf:
            fields = line.strip().split("\t")
            kraken_classified = bool(fields[idx["kraken_classified"]])
            kraken_assigned_host_virus = bool(fields[idx["kraken_assigned_host_virus"]])
            adj_score = float(fields[idx["bowtie2_length_normalized_score_max"]])
            # Discard if assigned to taxon other than host-infecting viruses
            if kraken_classified and not kraken_assigned_host_virus:
                continue
            # Discard if normalized Bowtie2 score is below threshold, unless assigned by Kraken to a host-infecting virus
            if !kraken_assigned_host_virus and adj_score < score_threshold:
                continue
            # Otherwise, write line to output
            outf.write("\t".join(fields) + "\n")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Filter putative virus reads based on Kraken and Bowtie2 results.")
    parser.add_argument("input_path", help="Path to input TSV file.")
    parser.add_argument("score_threshold", help="Length-normalized alignment score threshold for filtering.")
    parser.add_argument("output_file", help="Path to output TSV.")
    args = parser.parse_args()
    input_path = args.input_path
    score_threshold = float(args.score_threshold)
    out_path = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Input TSV file: {}".format(input_path))
    print_log("Score threshold: {}".format(score_threshold))
    print_log("Output TSV file: {}".format(out_path))
    # Run labeling function
    print_log("Filtering putative virus reads...")
    filter_virus_reads(input_path, score_threshold, out_path)
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
