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

def rename_columns(input_path, input_fields, output_fields, out_path):
    """Rename columns in TSV file."""
    if len(input_fields) != len(output_fields):
        raise ValueError("Input and output field lists must be the same length.")
    with open_by_suffix(input_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Read and handle header line
        headers_in = inf.readline().strip().split("\t")
        for field in input_fields:
            if field not in headers_in:
                raise ValueError(f"Input field not found in file header: {field}")
        headers_out = headers_in.copy()
        for i in range(len(input_fields)):
            headers_out[headers_in.index(input_fields[i])] = output_fields[i]
        header_line = "\t".join(headers_out)
        outf.write(header_line + "\n")
        # Write entire remainder of input file to output
        for line in inf:
            outf.write(line)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Rename one or more columns in a TSV file.")
    parser.add_argument("input_path", help="Path to input TSV file.")
    parser.add_argument("input_fields", help="Comma-separated list of input field names.")
    parser.add_argument("output_fields", help="Comma-separated list of output field names.")
    parser.add_argument("output_file", help="Path to output TSV.")
    args = parser.parse_args()
    input_path = args.input_path
    input_fields = args.input_fields.split(",")
    output_fields = args.output_fields.split(",")
    out_path = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Input TSV file: {}".format(input_path))
    print_log("Input fields: {}".format(input_fields))
    print_log("Output fields: {}".format(output_fields))
    print_log("Output TSV file: {}".format(out_path))
    # Run labeling function
    print_log("Renaming columns in TSV...")
    rename_columns(input_path, input_fields, output_fields, out_path)
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
