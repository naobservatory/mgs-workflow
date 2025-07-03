#!/usr/bin/env python3

"""
Add a prefix to column names in a TSV file.
Supports include/exclude modes for flexible column selection.
"""

# =======================================================================
# Preamble
# =======================================================================

import logging
import argparse
import csv
import sys
import gzip
import bz2
from typing import TextIO, List
from datetime import datetime, timezone


# Configure logging
class UTCFormatter(logging.Formatter):
    """
    Custom logging formatter that displays timestamps in UTC.

    Returns:
        Formatted log timestamps in UTC timezone
    """

    def formatTime(self, record, datefmt=None):
        """
        Format log timestamps in UTC timezone.

        Args:
            record: LogRecord object containing timestamp data
            datefmt: Optional date format string (unused)

        Returns:
            Formatted timestamp string in UTC
        """
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime("%Y-%m-%d %H:%M:%S UTC")


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter("[%(asctime)s] %(message)s")
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)

# =======================================================================
# I/O functions
# =======================================================================


def open_by_suffix(filename: str, mode: str = "r"):
    """
    Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files.

    Args:
        filename (str): Path to file to open
        mode (str): File open mode (default "r")

    Returns:
        File handle appropriate for the file compression type
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode + "t")
    elif filename.endswith(".bz2"):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


def parse_args() -> argparse.Namespace:
    """
    Parse and return command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Add prefix to column names in TSV file with include/exclude mode."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=lambda f: open_by_suffix(f, "r"),
        required=True,
        help="Input TSV file. Supports .gz and .bz2 compression.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=lambda f: open_by_suffix(f, "w"),
        required=True,
        help="Output TSV file. Supports .gz and .bz2 compression.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="Prefix to add to column names.",
    )
    parser.add_argument(
        "-c",
        "--columns",
        type=str,
        default="",
        help="Comma-separated list of column names. Empty string means no columns.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["include", "exclude"],
        required=True,
        help="Mode: 'include' to prefix only listed columns, 'exclude' to prefix all except listed columns.",
    )
    return parser.parse_args()


# =======================================================================
# Column functions
# =======================================================================


def parse_column_list(columns_str: str) -> List[str]:
    """
    Parse comma-separated column list, handling empty strings.

    Args:
        columns_str (str): Comma-separated list of column names

    Returns:
        List[str]: List of column names (empty list if input is empty)
    """
    if not columns_str or columns_str.isspace():
        return []
    return [col.strip() for col in columns_str.split(",") if col.strip()]


def apply_prefix_to_headers(
    headers: List[str], prefix: str, columns: List[str], mode: str, logger: logging.Logger
) -> List[str]:
    """
    Apply prefix to headers based on mode and column list.

    Args:
        headers (List[str]): Original column headers
        prefix (str): Prefix to add
        columns (List[str]): List of column names to include/exclude
        mode (str): 'include' or 'exclude'
        logger (logging.Logger): Logger instance

    Returns:
        List[str]: Headers with prefix applied according to rules
    """
    # Validate that specified columns exist in headers
    if columns:
        missing_columns = [col for col in columns if col not in headers]
        if missing_columns:
            error_msg = f"Specified columns not found in header: {missing_columns}. Available columns: {headers}"
            logger.error(error_msg)
            raise ValueError(error_msg)
    
    new_headers = []
    prefixed_count = 0
    
    for header in headers:
        if mode == "include":
            # Include mode: prefix only if column is in the list
            if header in columns:
                new_headers.append(f"{prefix}{header}")
                prefixed_count += 1
            else:
                new_headers.append(header)
        else:  # mode == "exclude"
            # Exclude mode: prefix only if column is NOT in the list
            if header not in columns:
                new_headers.append(f"{prefix}{header}")
                prefixed_count += 1
            else:
                new_headers.append(header)
    
    logger.info(f"Applied prefix '{prefix}' to {prefixed_count} out of {len(headers)} columns")
    return new_headers


# =======================================================================
# Main functions
# =======================================================================


def stream_and_prefix_tsv(
    input_file: TextIO,
    output_file: TextIO,
    prefix: str,
    columns: List[str],
    mode: str,
    logger: logging.Logger,
) -> None:
    """
    Stream TSV file and add prefix to column names based on mode.

    Args:
        input_file (TextIO): Input file handle
        output_file (TextIO): Output file handle
        prefix (str): Prefix to add to column names
        columns (List[str]): List of column names to include/exclude
        mode (str): 'include' or 'exclude'
        logger (logging.Logger): Logger instance
    """
    logger.info("Initializing streaming TSV reader and writer")
    reader = csv.reader(input_file, delimiter="\t")
    writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
    
    # Read and process header
    try:
        header = next(reader)
    except StopIteration:
        logger.warning("Input file is empty - writing empty output file")
        return
    
    logger.info(f"Input header has {len(header)} columns")
    logger.info(f"Mode: {mode}, Columns to {mode}: {columns if columns else 'none'}")
    
    # Apply prefix to headers
    new_header = apply_prefix_to_headers(header, prefix, columns, mode, logger)
    writer.writerow(new_header)
    logger.info("Header processed and written to output")
    
    # Stream remaining rows unchanged
    row_count = 0
    logger.info("Starting row-by-row streaming")
    
    for row_num, row in enumerate(reader, start=1):
        row_count += 1
        
        # Log progress for large files
        if row_count % 100000 == 0:
            logger.info(f"Processed {row_count} rows")
        
        # Validate row length
        if len(row) != len(header):
            logger.error(f"Row {row_num}: has {len(row)} columns but header has {len(header)} columns")
            raise ValueError(f"Row {row_num}: Column count mismatch")
        
        writer.writerow(row)
    
    # Final statistics
    logger.info(f"Streaming completed. Processed {row_count} data rows plus header")


def main() -> None:
    """Main function to execute the prefix addition process."""
    logger.info("Starting TSV column prefix addition")
    
    try:
        # Parse arguments
        args = parse_args()
        logger.info("Parsed command-line arguments successfully")
        
        # Parse column list
        columns = parse_column_list(args.columns)
        logger.info(f"Parsed {len(columns)} columns from input: {columns}")
        
        try:
            # Stream and process TSV
            logger.info("Starting streaming TSV processing...")
            stream_and_prefix_tsv(
                args.input,
                args.output,
                args.prefix,
                columns,
                args.mode,
                logger,
            )
            logger.info("Streaming processing completed successfully")
        finally:
            # Close file handles if they're not stdin/stdout
            if args.input != sys.stdin:
                args.input.close()
            if args.output != sys.stdout:
                args.output.close()
    
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        sys.exit(1)
    
    logger.info("TSV column prefix addition completed")


if __name__ == "__main__":
    main()
