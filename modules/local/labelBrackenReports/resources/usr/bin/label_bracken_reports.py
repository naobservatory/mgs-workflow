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

def label_bracken_report(input_path, sample_name, out_path):
    """Label Bracken report with header line and sample name."""
    with open_by_suffix(input_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Check header line has correct headers
        headers_exp = ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads",
                       "added_reads", "new_est_reads", "fraction_total_reads"]
        headers_obs = inf.readline().strip().split("\t")
        if headers_obs != headers_exp:
            print_log("Expected headers: {}".format(headers_exp))
            print_log("Observed headers: {}".format(headers_obs))
            raise ValueError("Header line does not match expected headers.")
        # Write header line to output
        headers_out = headers_exp + ["sample"]
        # Rename taxonomy_id to taxid and taxonomy_lvl to rank
        header_line = "\t".join(headers_out)
        outf.write(header_line + "\n")
        # Add sample name to each subsequent line and write to output
        for line in inf:
            outf.write(line.strip() + "\t" + sample_name + "\n")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Label Bracken reports with sample name.")
    parser.add_argument("bracken_report", help="Path to Bracken report file.")
    parser.add_argument("sample_name", help="Name of the sample.")
    parser.add_argument("output_file", help="Path to output TSV.")
    args = parser.parse_args()
    input_path = args.bracken_report
    sample_name = args.sample_name
    out_path = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Input Bracken report file: {}".format(input_path))
    print_log("Sample name: {}".format(sample_name))
    print_log("Output TSV file: {}".format(out_path))
    # Run labeling function
    print_log("Labeling Bracken report...")
    label_bracken_report(input_path, sample_name, out_path)
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
