#!/usr/bin/env python

"""
ABOUTME: Filters joined TSV from LCA and Bowtie2 data to keep only primary alignments.
ABOUTME: Removes rows where is_secondary=True, keeping only primary alignment data.
"""

import logging
import argparse
import csv
import sys
from typing import TextIO, List
from datetime import datetime, timezone


class UTCFormatter(logging.Formatter):
    """Custom formatter to display timestamps in UTC."""
    def formatTime(self, record, datefmt=None):
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime('%Y-%m-%d %H:%M:%S UTC')


def setup_logging() -> logging.Logger:
    """
    Set up logging configuration.
    
    Returns:
        logging.Logger: Configured logger instance.
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = UTCFormatter('[%(asctime)s] %(message)s')
    handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger


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
        "-i", "--input", 
        type=argparse.FileType('r'), 
        default=sys.stdin,
        help="Input joined TSV file (default: stdin)."
    )
    parser.add_argument(
        "-o", "--output", 
        type=argparse.FileType('w'), 
        default=sys.stdout,
        help="Output filtered TSV file (default: stdout)."
    )
    return parser.parse_args()


def validate_input_columns(header: List[str], logger: logging.Logger) -> int:
    """
    Validate that the input header contains required columns.
    
    Args:
        header: List of column names from CSV header.
        logger: Logger instance for logging messages.
        
    Returns:
        int: Index of the is_secondary column.
        
    Raises:
        ValueError: If required columns are missing.
    """
    required_columns = ['is_secondary']
    missing_columns = [col for col in required_columns if col not in header]
    
    if missing_columns:
        error_msg = f"Missing required columns: {missing_columns}"
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    is_secondary_index = header.index('is_secondary')
    logger.info(f"Input validation passed. Found {len(header)} columns.")
    return is_secondary_index


def convert_to_boolean(value: str) -> bool:
    """
    Convert string value to boolean.
    
    Args:
        value: String value to convert.
        
    Returns:
        bool: Converted boolean value.
        
    Raises:
        ValueError: If value cannot be converted to boolean.
    """
    if value.lower() in ('true', '1'):
        return True
    elif value.lower() in ('false', '0'):
        return False
    else:
        raise ValueError(f"Cannot convert '{value}' to boolean")


def filter_primary_alignments(input_file: TextIO, output_file: TextIO, logger: logging.Logger) -> None:
    """
    Filter TSV file to keep only primary alignments.
    
    Args:
        input_file: Input file handle.
        output_file: Output file handle.
        logger: Logger instance for logging messages.
    """
    reader = csv.reader(input_file, delimiter='\t')
    writer = csv.writer(output_file, delimiter='\t')
    
    # Read and validate header
    header = next(reader)
    is_secondary_index = validate_input_columns(header, logger)
    writer.writerow(header)
    
    initial_rows = 0
    filtered_rows = 0
    conversion_errors = 0
    
    for row in reader:
        initial_rows += 1
        
        try:
            is_secondary_value = row[is_secondary_index]
            is_secondary = convert_to_boolean(is_secondary_value)
            
            # Keep only primary alignments (is_secondary == False)
            if not is_secondary:
                writer.writerow(row)
                filtered_rows += 1
                
        except (ValueError, IndexError) as e:
            conversion_errors += 1
            logger.warning(f"Row {initial_rows}: {str(e)}")
            continue
    
    removed_rows = initial_rows - filtered_rows
    
    logger.info(f"Input contains {initial_rows} rows.")
    logger.info(f"Filtered to {filtered_rows} primary alignments.")
    logger.info(f"Removed {removed_rows} secondary alignments.")
    
    if conversion_errors > 0:
        logger.warning(f"Encountered {conversion_errors} rows with conversion errors.")


def main() -> None:
    """Main function to execute the filtering process."""
    logger = setup_logging()
    logger.info("Starting primary alignment filtering.")
    
    try:
        # Parse arguments
        args = parse_args()
        logger.info("Parsed command-line arguments successfully.")
        
        # Filter primary alignments
        logger.info("Processing TSV data...")
        filter_primary_alignments(args.input, args.output, logger)
        logger.info("Output written successfully.")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        sys.exit(1)
    finally:
        # Close file handles if they're not stdin/stdout
        if hasattr(args, 'input') and args.input != sys.stdin:
            args.input.close()
        if hasattr(args, 'output') and args.output != sys.stdout:
            args.output.close()
        
        logger.info("Primary alignment filtering completed.")


if __name__ == "__main__":
    main()
