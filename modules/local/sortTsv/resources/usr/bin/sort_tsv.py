#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import subprocess
import logging
from datetime import datetime
import time

# Set up logging with UTC time
class UTCFormatter(logging.Formatter):
    converter = time.gmtime

# Configure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = UTCFormatter('[%(asctime)s UTC] %(levelname)s: %(message)s', 
                         datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)

def open_by_suffix(filename, mode="r"):
    """Open a file with appropriate mode based on suffix."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)

def sort_tsv_file(input_file, output_file, sort_field):
    """Sort a TSV file by a specific column header."""
    # Create temporary file for data rows
    temp_sorted = "temp_sorted.txt"
    
    # Process input file line by line
    with open_by_suffix(input_file) as f:
        # Read the header first
        header = f.readline().strip()
        
        # Check if file is empty (no header)
        if not header:
            logger.warning(f"Input file {input_file} is empty. Creating empty output file.")
            # Create empty output file by writing an empty string
            with open_by_suffix(output_file, "w") as out:
                out.write("")
            return
            
        # Split header to find sort column
        header_fields = header.split('\t')
        
        try:
            # Find the index of the sort field (0-based)
            col_index = header_fields.index(sort_field)
        except ValueError:
            msg = f"Could not find sort field in input header: '{sort_field}', {header_fields}"
            logger.error(msg)
            raise ValueError(msg)
            
        logger.info(f"Found '{sort_field}' in column {col_index+1}. Sorting by column {col_index+1}.")
        
        # Check if there are any data rows by reading the first data line
        first_data_line = f.readline()
        if not first_data_line:
            logger.warning("No data rows to sort, writing header only")
            # Just write the header to output
            with open_by_suffix(output_file, "w") as out:
                out.write(header + "\n")
            return
        
        # We have at least one data line, write it and all subsequent lines to temp file
        with open(temp_sorted, "w") as temp:
            temp.write(first_data_line)  # Write the first data line we already read
            # Write remaining lines
            for line in f:
                temp.write(line)
    
    # Sort the data using the sort command
    try:
        sort_cmd = ["sort", "-t", "\t", "-k", f"{col_index+1},{col_index+1}", temp_sorted]
        sorted_data = subprocess.run(sort_cmd, capture_output=True, text=True, check=True).stdout
        
        # Write header and sorted data to output
        with open_by_suffix(output_file, "w") as out:
            out.write(header + "\n")
            out.write(sorted_data)
    except subprocess.CalledProcessError as e:
        logger.error(f"Sort command failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Error writing sorted data: {e}")
        raise
    finally:
        # Clean up temporary file
        if os.path.exists(temp_sorted):
            os.remove(temp_sorted)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Sort a TSV file by a specific column header.")
    parser.add_argument("input_file", help="Input TSV file path")
    parser.add_argument("sort_field", help="Column header to sort by")
    parser.add_argument("output_file", help="Output TSV file path")
    
    try:
        args = parser.parse_args()
        
        # Start time tracking
        start_time = datetime.now()
        logger.info("Starting sort process")
        
        # Log parameters
        logger.info(f"Input file: {args.input_file}")
        logger.info(f"Sort field: {args.sort_field}")
        logger.info(f"Output file: {args.output_file}")
        
        # Sort the TSV file
        sort_tsv_file(args.input_file, args.output_file, args.sort_field)
        
        # Calculate elapsed time
        elapsed_time = (datetime.now() - start_time).total_seconds()
        logger.info(f"Completed in {elapsed_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Failed to sort TSV file: {e}")
        raise

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.exit(1)
