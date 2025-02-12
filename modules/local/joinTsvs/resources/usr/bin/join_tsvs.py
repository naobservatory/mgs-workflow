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

def write_line(line_list, output_file):
    """Write a line to the output file."""
    output_file.write("\t".join(line_list) + "\n")

def fill_right(row_1, placeholder_2):
    """Fill in fields from file 2 with placeholder values."""
    return row_1 + placeholder_2

def fill_left(placeholder_1, row_2, id_2, field_index_1, field_index_2):
    """Fill in fields from file 1 with placeholder values."""
    # Put join field in appropriate position in columns from file 1
    if field_index_1 == 0:
        merged_row = [id_2] + placeholder_1
    else:
        merged_row = placeholder_1[:field_index_1] + [id_2] + placeholder_1[field_index_1:]
    # Add remainder of row_2 excluding join field
    merged_row += [val for i, val in enumerate(row_2) if i != field_index_2]
    return merged_row

def join_tsvs(input_path_1, input_path_2, field, join_type, output_file):
    """Join two TSV files linewise on a shared column."""
    # Open files
    with open_by_suffix(input_path_1, "r") as file_1, open_by_suffix(input_path_2, "r") as file_2, open_by_suffix(output_file, "w") as output:
        # Read header lines
        header_1 = file_1.readline().strip().split("\t")
        header_2 = file_2.readline().strip().split("\t")
        # Verify that the field exists in both files
        if field not in header_1:
            raise ValueError(f"Join field missing from file 1 ('{field}', {input_path_1}).")
        if field not in header_2:
            raise ValueError(f"Join field missing from file 2 ('{field}', {input_path_2}).")
        # Get index of join field
        field_index_1 = header_1.index(field)
        field_index_2 = header_2.index(field)
        # Write new header
        merged_header = header_1 + [h for i, h in enumerate(header_2) if i != field_index_2]
        write_line(merged_header, output)
        # Define placeholders for non-inner joins
        placeholder_file1 = ["NA"] * (len(header_1) - 1)
        placeholder_file2 = ["NA"] * (len(header_2) - 1)
        # Initialise trackers for verifying files are sorted
        last_id_file1 = None
        last_id_file2 = None
        # Get first line from each file
        line_1 = file_1.readline()
        line_2 = file_2.readline()
        # Iterate until we exhaust either file
        while line_1 and line_2:
            # Get IDs and verify sorting
            row_1 = line_1.rstrip("\n").split("\t")
            row_2 = line_2.rstrip("\n").split("\t")
            id_1 = row_1[field_index_1]
            id_2 = row_2[field_index_2]
            #print(id_1, id_2, last_id_file1, last_id_file2)
            if last_id_file1 is not None and id_1 < last_id_file1:
                raise ValueError(f"File 1 is not sorted: encountered ID {id_1} after {last_id_file1} ({input_path_1}).")
            if last_id_file2 is not None and id_2 < last_id_file2:
                raise ValueError(f"File 2 is not sorted: encountered ID {id_2} after {last_id_file2} ({input_path_2}).")
            last_id_file1 = id_1
            last_id_file2 = id_2
            if id_1 == id_2:
                # Regardless of join type, if IDs match, write to output and advance both files
                merged_row = row_1 + [val for i, val in enumerate(row_2) if i != field_index_2]
                write_line(merged_row, output)
                line_1 = file_1.readline()
                line_2 = file_2.readline()
            elif id_1 < id_2:
                # File 2 is missing ID present in file 1
                if join_type in ("left", "outer"):
                    merged_row = fill_right(row_1, placeholder_file2)
                    write_line(merged_row, output)
                line_1 = file_1.readline()
            else:
                # File 1 is missing ID present in file 2
                if join_type in ("right", "outer"):
                    merged_row = fill_left(placeholder_file1, row_2, id_2, field_index_1, field_index_2)
                    write_line(merged_row, output)
                line_2 = file_2.readline()
        # Read out file 1
        while line_1:
            row_1 = line_1.rstrip("\n").split("\t")
            id_1 = row_1[field_index_1]
            if last_id_file1 is not None and id_1 < last_id_file1:
                raise ValueError(f"File 1 is not sorted: encountered ID {id_1} after {last_id_file1} ({input_path_1}).")
            last_id_file1 = id_1
            if join_type in ("left", "outer"):
                merged_row = fill_right(row_1, placeholder_file2)
                write_line(merged_row, output)
            line_1 = file_1.readline()
        # Read out file 2
        while line_2:
            row_2 = line_2.rstrip("\n").split("\t")
            id_2 = row_2[field_index_2]
            if last_id_file2 is not None and id_2 < last_id_file2:
                raise ValueError(f"File 2 is not sorted: encountered ID {id_2} after {last_id_file2} ({input_path_2}).")
            last_id_file2 = id_2
            if join_type in ("right", "outer"):
                merged_row = fill_left(placeholder_file1, row_2, id_2, field_index_1, field_index_2)
                write_line(merged_row, output)
            line_2 = file_2.readline()

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Join two TSV files linewise on a shared column.")
    parser.add_argument("tsv1", help="Path to first TSV file.")
    parser.add_argument("tsv2", help="Path to second TSV file.")
    parser.add_argument("field", type=str, help="Field to join on.")
    parser.add_argument("join_type", type=str, help="Type of join to perform (inner, left, right, outer).")
    parser.add_argument("output_file", type=str, help="Path to output TSV file.")
    args = parser.parse_args()
    # Assign variables
    input_path_1 = args.tsv1
    input_path_2 = args.tsv2
    field = args.field
    join_type = args.join_type
    output_file = args.output_file
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Input TSV file 1: {}".format(input_path_1))
    print_log("Input TSV file 2: {}".format(input_path_2))
    print_log("Field to join on: {}".format(field))
    print_log("Join type: {}".format(join_type))
    print_log("Output file: {}".format(output_file))
    # Run joining function
    join_tsvs(input_path_1, input_path_2, field, join_type, output_file)
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
