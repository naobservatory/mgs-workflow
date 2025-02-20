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

def process_headers(header_1, header_2, field):
    """Read, join and check headers from both files."""
    # Verify that the field exists in both files
    if field not in header_1:
        raise ValueError(f"Join field missing from file 1 ('{field}').")
    if field not in header_2:
        raise ValueError(f"Join field missing from file 2 ('{field}').")
    # Check for duplicate field names across files
    for h in header_1:
        if h != field and h in header_2:
            msg = f"Duplicate field name found across both files: '{h}'."
            raise ValueError(msg)
    for h in header_2:
        if h != field and h in header_1:
            msg = f"Duplicate field name found across both files: '{h}'."
            raise ValueError(msg)
    # Get index of join field
    field_index_1 = header_1.index(field)
    field_index_2 = header_2.index(field)
    # Prepare new header
    merged_header = header_1 + \
            [h for i, h in enumerate(header_2) if i != field_index_2]
    return merged_header, field_index_1, field_index_2

def get_line_id(file, field_index):
    """Get the next line ID from a file."""
    line = file.readline()
    row = line.rstrip("\n").split("\t") if line else None
    id = row[field_index] if row else None
    return line, row, id

def check_sorting(id_curr, id_next, file_str, input_path):
    """Check that the file is sorted and raise an error if not."""
    if id_next is not None and id_curr > id_next:
        msg = f"File {file_str} is not sorted: encountered ID {id_curr} after {id_next} ({input_path})."
        raise ValueError(msg)

def join_tsvs(input_path_1, input_path_2, field, join_type, output_file):
    """Join two TSV files linewise on a shared column."""
    # Open files
    with open_by_suffix(input_path_1, "r") as file_1, open_by_suffix(input_path_2, "r") as file_2, open_by_suffix(output_file, "w") as output:
        # Read header lines
        header_1 = file_1.readline().strip().split("\t")
        header_2 = file_2.readline().strip().split("\t")
        # Process and write header
        merged_header, field_index_1, field_index_2 = \
                process_headers(header_1, header_2, field)
        write_line(merged_header, output)
        # Define placeholders for non-inner joins
        placeholder_file1 = ["NA"] * (len(header_1) - 1)
        placeholder_file2 = ["NA"] * (len(header_2) - 1)
        # Get first two lines from each file
        line_1_curr, row_1_curr, id_1_curr = get_line_id(file_1, field_index_1)
        line_1_next, row_1_next, id_1_next = get_line_id(file_1, field_index_1)
        line_2_curr, row_2_curr, id_2_curr = get_line_id(file_2, field_index_2)
        line_2_next, row_2_next, id_2_next = get_line_id(file_2, field_index_2)
        # Iterate until we exhaust either file
        while line_1_curr and line_2_curr:
            # Verify that files are sorted
            check_sorting(id_1_curr, id_1_next, "1", input_path_1)
            check_sorting(id_2_curr, id_2_next, "2", input_path_2)
            # If IDs match, write to output and check for one-to-many join
            if id_1_curr == id_2_curr:
                merged_row = row_1_curr + [val for i, val in enumerate(row_2_curr) if i != field_index_2]
                write_line(merged_row, output)
                if id_1_curr != id_1_next and id_2_curr != id_2_next:
                    line_1_curr, row_1_curr, id_1_curr = line_1_next, row_1_next, id_1_next
                    line_1_next, row_1_next, id_1_next = get_line_id(file_1, field_index_1)
                    line_2_curr, row_2_curr, id_2_curr = line_2_next, row_2_next, id_2_next
                    line_2_next, row_2_next, id_2_next = get_line_id(file_2, field_index_2)
                elif id_1_curr != id_1_next:
                    line_2_curr, row_2_curr, id_2_curr = line_2_next, row_2_next, id_2_next
                    line_2_next, row_2_next, id_2_next = get_line_id(file_2, field_index_2)
                elif id_2_curr != id_2_next:
                    line_1_curr, row_1_curr, id_1_curr = line_1_next, row_1_next, id_1_next
                    line_1_next, row_1_next, id_1_next = get_line_id(file_1, field_index_1)
                else:
                    ids = f"{id_1_curr}, {id_1_next}, {id_2_curr}, {id_2_next}"
                    msg = f"Unsupported many-to-many join detected for ID {id_1_curr} ({ids})."
                    raise ValueError(msg)
            # If IDs don't match, handle on the basis of join type
            elif id_1_curr < id_2_curr:
                # File 2 is missing ID present in file 1
                if join_type == "strict":
                    msg = f"Strict join failed: ID {id_1_curr} missing from file 2."
                    raise ValueError(msg)
                elif join_type in ("left", "outer"):
                    merged_row = fill_right(row_1_curr, placeholder_file2)
                    write_line(merged_row, output)
                line_1_curr, row_1_curr, id_1_curr = line_1_next, row_1_next, id_1_next
                line_1_next, row_1_next, id_1_next = get_line_id(file_1, field_index_1)
            else:
                # File 1 is missing ID present in file 2
                if join_type == "strict":
                    msg = f"Strict join failed: ID {id_2_curr} missing from file 1."
                    raise ValueError(msg)
                if join_type in ("right", "outer"):
                    merged_row = fill_left(placeholder_file1, row_2_curr, id_2_curr, field_index_1, field_index_2)
                    write_line(merged_row, output)
                line_2_curr, row_2_curr, id_2_curr = line_2_next, row_2_next, id_2_next
                line_2_next, row_2_next, id_2_next = get_line_id(file_2, field_index_2)
        # Read out file 1
        while line_1_curr:
            check_sorting(id_1_curr, id_1_next, "1", input_path_1)
            if join_type == "strict":
                msg = f"Strict join failed: ID {id_1_curr} missing from file 2."
                raise ValueError(msg)
            if join_type in ("left", "outer"):
                merged_row = fill_right(row_1_curr, placeholder_file2)
                write_line(merged_row, output)
            line_1_curr, row_1_curr, id_1_curr = line_1_next, row_1_next, id_1_next
            line_1_next, row_1_next, id_1_next = get_line_id(file_1, field_index_1)
        # Read out file 2
        while line_2_curr:
            check_sorting(id_2_curr, id_2_next, "2", input_path_2)
            if join_type == "strict":
                msg = f"Strict join failed: ID {id_2_curr} missing from file 1."
                raise ValueError(msg)
            if join_type in ("right", "outer"):
                merged_row = fill_left(placeholder_file1, row_2_curr, id_2_curr, field_index_1, field_index_2)
                write_line(merged_row, output)
            line_2_curr, row_2_curr, id_2_curr = line_2_next, row_2_next, id_2_next
            line_2_next, row_2_next, id_2_next = get_line_id(file_2, field_index_2)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Join two TSV files linewise on a shared column.")
    parser.add_argument("tsv1", help="Path to first TSV file.")
    parser.add_argument("tsv2", help="Path to second TSV file.")
    parser.add_argument("field", type=str, help="Field to join on.")
    parser.add_argument("join_type", type=str, help="Type of join to perform (inner, left, right, outer, strict).")
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
    # Verify that join type is valid
    if join_type not in ("inner", "left", "right", "outer", "strict"):
        raise ValueError("Invalid join type. Must be one of 'inner', 'left', 'right', 'outer', or 'strict'.")
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
