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

def add_sample_column(input_path, sample_name, sample_column, out_path):
    """Add sample name column to TSV file."""
    with open_by_suffix(input_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Read and handle header line
        headers_in = inf.readline().strip().split("\t")
        if sample_column in headers_in:
            raise ValueError(f"Sample name column already exists: {sample_column}")
        headers_out = headers_in + [sample_column]
        header_line = "\t".join(headers_out)
        outf.write(header_line + "\n")
        # Add sample name to each subsequent line and write to output
        for line in inf:
            outf.write(line.strip() + "\t" + sample_name + "\n")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Add sample name column to TSV file.")
    parser.add_argument("input_path", help="Path to input TSV file.")
    parser.add_argument("sample_name", help="Name of the sample.")
    parser.add_argument("sample_column", help="Name of the sample name column.")
    parser.add_argument("output_file", help="Path to output TSV.")
    args = parser.parse_args()
    input_path = args.bracken_report
    sample_name = args.sample_name
    sample_column = args.sample_column
    out_path = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Input Bracken report file: {}".format(input_path))
    print_log("Sample name: {}".format(sample_name))
    print_log("Sample name column: {}".format(sample_column))
    print_log("Output TSV file: {}".format(out_path))
    # Run labeling function
    print_log("Adding sample name column to TSV...")
    add_sample_column(input_path, sample_name, sample_column, out_path)
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
