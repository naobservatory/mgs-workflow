#!/usr/bin/env python

"""
Given a TSV column with an initial header line,
return a new TSV containing only the specified columns.
"""

#=======================================================================
# Import modules
#=======================================================================

import logging
from datetime import datetime, timezone
import argparse
import time
import gzip
import bz2

#=======================================================================
# Configure logging
#=======================================================================

class UTCFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime('%Y-%m-%d %H:%M:%S UTC')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter('[%(asctime)s] %(message)s')
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)

#=======================================================================
# I/O functions
#=======================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    # Create parser
    desc = "Given a TSV column with an initial header line, " \
           "return a new TSV containing only the specified columns."
    parser = argparse.ArgumentParser(description=desc)
    # Add arguments
    parser.add_argument("--input", "-i", help="Path to input TSV.")
    parser.add_argument("--output", "-o", help="Path to output TSV.")
    parser.add_argument("--fields", "-f", help="Comma-separated list of fields to select.")
    # Return parsed arguments
    return parser.parse_args()

def open_by_suffix(filename, mode="r"):
    logger.debug(f"\tOpening file object: {filename}")
    logger.debug(f"\tOpening mode: {mode}")
    logger.debug(f"\tGZIP mode: {filename.endswith('.gz')}")
    logger.debug(f"\tBZ2 mode: {filename.endswith('.bz2')}")
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

#=======================================================================
# TSV processing functions
#=======================================================================

def get_header_index(headers: list[str], field: str) -> int:
    """Get the index of a field in a header line."""
    try:
        return headers.index(field)
    except ValueError:
        raise ValueError(f"Field not found in header: {field}")
    
def join_line(inputs: list[str]) -> str:
    """Join a list of strings with tabs followed by a newline."""
    return "\t".join(inputs) + "\n"

def subset_line(inputs: list[str], indices: list[int]) -> list[str]:
    """Subset a list of strings to a list of indices."""
    return [inputs[i] for i in indices]

def select_columns(input_path: str, output_path: str, fields: list[str]) -> None:
    """Select columns in TSV file."""
    with open_by_suffix(input_path) as inf, open_by_suffix(output_path, "w") as outf:
        # Read the header line
        header_line = inf.readline().strip()
        if not header_line:
            # No header = selected fields are missing
            raise ValueError("No header to select fields from.")
        headers_in = header_line.split("\t")
        # Get indices of selected fields
        indices = [get_header_index(headers_in, field) for field in fields]
        # Write new header
        new_header_line = join_line(subset_line(headers_in, indices))
        outf.write(new_header_line)
        # Subset and write body lines of input file
        for line in inf:
            line = line.strip()
            if not line:
                continue
            line_split_subset = subset_line(line.split("\t"), indices)
            outf.write(join_line(line_split_subset))

#=======================================================================
# Main function
#=======================================================================

def main():
    logger.info("Initializing script.")
    start_time = time.time()
    # Parse arguments
    logger.info("Parsing arguments.")
    args = parse_args()
    input_path = args.input
    output_path = args.output
    logger.info(f"Input TSV file: {input_path}")
    logger.info(f"Output TSV file: {output_path}")
    # Extract fields to select from header
    fields = args.fields.split(",")
    logger.info(f"Selected fields: {fields}")
    if len(fields) == 0:
        raise ValueError("No fields to select.")
    # Select columns from input TSV
    logger.info("Selecting columns from input TSV.")
    select_columns(input_path, output_path, fields)
    # Finish time tracking
    logger.info("Script completed successfully.")
    end_time = time.time()
    logger.info(f"Total time elapsed: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()