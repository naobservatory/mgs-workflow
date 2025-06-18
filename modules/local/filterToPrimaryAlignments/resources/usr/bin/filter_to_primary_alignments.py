#!/usr/bin/env python3

"""
Filters joined TSV from LCA and Bowtie2 data to keep only primary alignments.
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
from typing import TextIO
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
        description="Filter joined TSV to keep only primary alignments (is_secondary=False)."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=lambda f: open_by_suffix(f, "r"),
        default=sys.stdin,
        help="Input joined TSV file (default: stdin). Supports .gz and .bz2 compression.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=lambda f: open_by_suffix(f, "w"),
        default=sys.stdout,
        help="Output filtered TSV file (default: stdout). Supports .gz and .bz2 compression.",
    )
    return parser.parse_args()


# =======================================================================
# Column functions
# =======================================================================


def validate_input_columns(header: list[str], logger: logging.Logger) -> int:
    """
    Validate that the input header contains required columns.

    Args:
        header (list[str]): list of column names from CSV header.
        logger (logging.Logger): Logger instance for logging messages.

    Returns:
        int: Index of the is_secondary column.

    Raises:
        ValueError: If required columns are missing.
    """
    required_columns = ["is_secondary"]
    missing_columns = [col for col in required_columns if col not in header]

    if missing_columns:
        error_msg = f"Missing required columns: {missing_columns}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    is_secondary_index = header.index("is_secondary")
    logger.info(f"Input validation passed. Found {len(header)} columns.")
    return is_secondary_index




# =======================================================================
# Main functions
# =======================================================================


def stream_and_filter_tsv(
    input_file: TextIO, output_file: TextIO, logger: logging.Logger
) -> None:
    """
    Stream TSV file and filter to keep only primary alignments using memory-efficient processing.
    
    Processes the input file row by row without loading the entire dataset into memory,
    making it suitable for large files.

    Args:
        input_file (TextIO): Input file handle.
        output_file (TextIO): Output file handle.
        logger (logging.Logger): Logger instance for logging messages.
    """
    logger.info("Initializing streaming TSV reader and writer")
    reader = csv.reader(input_file, delimiter="\t")
    writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")

    # Read and validate header
    try:
        header = next(reader)
    except StopIteration:
        logger.error("Input file is empty or contains no header")
        raise ValueError("Input file is empty or contains no header")
    
    is_secondary_index = validate_input_columns(header, logger)
    writer.writerow(header)
    logger.info("Header processed and written to output")

    # Initialize counters
    initial_rows = 0
    filtered_rows = 0

    logger.info("Starting row-by-row streaming processing")
    
    # Stream and process each row individually
    for row_num, row in enumerate(reader, start=1):
        initial_rows += 1
        
        # Log progress for large files
        if initial_rows % 100000 == 0:
            logger.info(f"Processed {initial_rows} rows, kept {filtered_rows} primary alignments")

        try:
            # Check row length to prevent IndexError
            if len(row) <= is_secondary_index:
                error_msg = f"Row {row_num}: Row has {len(row)} columns but expected at least {is_secondary_index + 1}"
                logger.error(error_msg)
                raise ValueError(error_msg)
                
            is_secondary_value = row[is_secondary_index]
            # Convert string representation to boolean if needed
            if isinstance(is_secondary_value, str):
                is_secondary = is_secondary_value.lower() in ('true', '1')
            else:
                is_secondary = bool(is_secondary_value)

            # Keep only primary alignments (is_secondary == False)
            if not is_secondary:
                writer.writerow(row)
                filtered_rows += 1

        except (ValueError, IndexError) as e:
            logger.error(f"Row {row_num}: {str(e)}")
            logger.error("Exiting due to data format error")
            raise

    # Final statistics
    removed_rows = initial_rows - filtered_rows

    logger.info("Streaming processing completed:")
    logger.info(f"  - Input rows processed: {initial_rows}")
    logger.info(f"  - Primary alignments kept: {filtered_rows}")
    logger.info(f"  - Secondary alignments removed: {removed_rows}")


def main() -> None:
    """Main function to execute the filtering process."""
    logger.info("Starting primary alignment filtering.")

    try:
        # Parse arguments
        args = parse_args()
        logger.info("Parsed command-line arguments successfully.")

        try:
            # Stream and filter primary alignments
            logger.info("Starting streaming TSV processing...")
            stream_and_filter_tsv(args.input, args.output, logger)
            logger.info("Streaming processing completed successfully.")
        finally:
            # Close file handles if they're not stdin/stdout
            if args.input != sys.stdin:
                args.input.close()
            if args.output != sys.stdout:
                args.output.close()

    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        sys.exit(1)

    logger.info("Primary alignment filtering completed.")


if __name__ == "__main__":
    main()
