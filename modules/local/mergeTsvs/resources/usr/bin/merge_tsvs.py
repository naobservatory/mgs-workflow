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

def label_kraken_report(input_path, sample_name, out_path):
    """Label Kraken report with header line and sample name."""
    with open_by_suffix(input_path) as inf, open_by_suffix(out_path, "w") as outf:
        # Add header line to output file
        headers = ["pc_reads_total", "n_reads_clade", "n_reads_direct",
                   "n_minimizers_total", "n_minimizers_distinct", "rank",
                   "taxid", "name", "sample"]
        header_line = "\t".join(headers)
        outf.write(header_line + "\n")
        # Add sample name to each line in input and write to output
        for line in inf:
            outf.write(line.strip() + "\t" + sample_name + "\n")

def read_header(infile) -> List[str]:
    """Read header from TSV file and return as list."""
    header_line = infile.readline()
    if not header_line:
        raise ValueError("File is empty.")
    header = header_line.strip().split("\t")
    return header

def map_headers(header: List[str], reference_header: List[str]) -> List[int]:
    """Generate mapping of columns in header to reference header."""
    header_mapping = {col: i for i, col in enumerate(header)}
    reference_mapping = [header_mapping[col] for col in reference_header]
    return reference_mapping

def check_headers(header: List[str], reference_header: List[str]) -> bool:
    """Check if fields in two headers match and raise error if not."""
    hset = set(header)
    rset = set(reference_header)
    header_difference = hset ^ rset
    if not header_difference:
        return True
    msg = f"Headers do not match:"
    if rset - hset:
        msg += f"\n\tMissing fields: {rset - hset}"
    if hset - rset:
        msg += f"\n\tExtra fields: {hset - rset}"
    raise ValueError(msg)

def concatenate_tsvs(input_files: List[str], out_path: str):
    """Concatenate multiple TSV files with matching headers."""
    with open_by_suffix(out_path, "w") as outf:
        # Parse first file and get header mapping
        with open_by_suffix(input_files[0]) as inf:
            print_log(f"\tProcessing first file: {input_files[0]}")
            # Extract header information
            reference_header = read_header(inf)
            # Write header to output file
            outf.write("\t".join(reference_header) + "\n")
            # Write contents of first file to output unchanged
            for line in inf:
                outf.write(line)
        # Parse remaining files and write to output
        print_log(f"\tProcessing other files:")
        for input_path in input_files[1:]:
            print_log(f"\t\tProcessing file: {input_path}")
            with open_by_suffix(input_path) as inf:
                # Extract header information
                header = read_header(inf)
                # Verify header fields match reference header
                check_headers(header, reference_header)
                # Generate mapping of header fields to reference header
                header_mapping = map_headers(header, reference_header)
                # Write contents of file to output with mapped fields
                for line in inf:
                    fields = line.strip().split("\t")
                    mapped_fields = [fields[i] for i in header_mapping]
                    outf.write("\t".join(mapped_fields) + "\n")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Concatenate multiple TSV files with matching headers.")
    parser.add_argument("input_files", nargs="+", help="Paths to input TSV files.",
                        metavar="INPUT", type=str)
    parser.add_argument("-o", "--output_file", help="Path to output TSV.",
                        type=str, required=True)
    args = parser.parse_args()
    input_files = args.input_files
    out_path = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Output TSV file: {}".format(out_path))
    print_log("Output TSV files:")
    for input_path in input_files:
        print_log(f"\t{input_path}")
    # Concatenate input files
    print_log("Concatenating TSV files...")
    concatenate_tsvs(input_files, out_path)
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
