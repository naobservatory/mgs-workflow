#!/usr/bin/env python

"""
Validate grouping file against virus hits file using directional validation.

This script implements asymmetric validation:
1. ENFORCE: Every sample in virus_hits MUST have corresponding grouping (data integrity)
2. ALLOW: Samples in grouping without virus hits (zero viral hits - filter out)

Usage:
    validate_grouping.py virus_hits.tsv grouping.tsv join_field validated_grouping.tsv.gz samples_without_vv_hits.tsv
"""

#=======================================================================
# Preamble
#=======================================================================

# Import libraries
import argparse
import logging
import gzip
import bz2
from collections import defaultdict
from datetime import datetime, timezone
from typing import TextIO

# Configure logging
class UTCFormatter(logging.Formatter):
    """Custom formatter that outputs timestamps in UTC format.
    
    This formatter ensures all log timestamps are in UTC for consistency
    across different timezones and environments.
    """
    
    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        """Format the time in UTC.
        
        Args:
            record: LogRecord instance containing the event data
            datefmt: Date format string (currently unused, UTC format is hardcoded)
        
        Returns:
            str: Formatted UTC timestamp string
        """
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

def open_by_suffix(filename: str, mode: str = "r") -> TextIO:
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files.
    
    Args:
        filename: Path to the file to open
        mode: File open mode (default: "r")
    
    Returns:
        TextIO: File handle for reading or writing
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.
    
    Args:
        None
    
    Returns:
        argparse.Namespace: Parsed command-line arguments containing:
            - virus_hits_file: Path to virus hits TSV file
            - grouping_file: Path to grouping TSV file
            - join_field: Field to validate on
            - validated_output: Path to validated grouping output file
            - samples_without_vv_hits: Path to samples without viral hits log file
    """
    parser = argparse.ArgumentParser(description="Validate grouping file against virus hits file")
    parser.add_argument("virus_hits_file", help="Path to virus hits TSV file")
    parser.add_argument("grouping_file", help="Path to grouping TSV file")
    parser.add_argument("join_field", help="Field to validate on (e.g., 'sample')")
    parser.add_argument("validated_output", help="Path to validated grouping output file (.tsv.gz)")
    parser.add_argument("samples_without_vv_hits", help="Path to samples without viral hits log file (.tsv)")
    return parser.parse_args()

#=======================================================================
# File reading functions
#=======================================================================

def read_sample_ids_from_file(file_path: str, join_field: str) -> tuple[set[str], list[str]]:
    """
    Read sample IDs from a TSV file.
    
    Args:
        file_path: Path to the TSV file to read
        join_field: Name of the field to extract sample IDs from
    
    Returns:
        tuple[set[str], list[str]]: A tuple containing:
            - Set of unique sample IDs found in the file
            - List of header field names from the file
    
    Raises:
        ValueError: If the join_field is not found in the file header
    """
    sample_ids = set()
    
    with open_by_suffix(file_path) as file_handle:
        header_line = file_handle.readline().strip()
        if not header_line:
            logger.warning(f"File {file_path} is empty")
            return set(), []
        header = header_line.split('\t')
        if join_field not in header:
            raise ValueError(f"Join field '{join_field}' not found in {file_path}. Available fields: {header}")
        join_field_index = header.index(join_field)
        logger.info(f"Reading sample IDs from {file_path}, join field '{join_field}' at index {join_field_index}")
        for line_num, line in enumerate(file_handle, start=2):
            line = line.strip()
            fields = line.split('\t')
            sample_id = fields[join_field_index]
            sample_ids.add(sample_id)
        logger.info(f"Found {len(sample_ids)} unique sample IDs in {file_path}")
        return sample_ids, header

#=======================================================================
# Validation functions
#=======================================================================

def validate_directional_integrity(virus_hits_ids: set[str], grouping_ids: set[str]) -> tuple[set[str], set[str]]:
    """
    Perform directional validation between virus hits and grouping files.
    
    This function implements asymmetric validation where samples with virus hits
    must have grouping (data integrity), but samples with grouping may not have
    virus hits (zero viral hits).
    
    Args:
        virus_hits_ids: Set of sample IDs from virus hits file
        grouping_ids: Set of sample IDs from grouping file
    
    Returns:
        tuple[set[str], set[str]]: A tuple containing:
            - missing_grouping_ids: Samples in virus_hits but not in grouping (ERROR)
            - zero_vv_ids: Samples in grouping but not in virus_hits (OK - zero viral hits)
    """
    # Find samples with virus hits but missing grouping (DATA INTEGRITY ERROR)
    missing_grouping = virus_hits_ids - grouping_ids
    # Find samples with grouping but no virus hits (ZERO VIRAL HITS - OK)
    zero_vv = grouping_ids - virus_hits_ids
    logger.info(f"Validation results:")
    logger.info(f"  - Samples with virus hits: {len(virus_hits_ids)}")
    logger.info(f"  - Samples with grouping: {len(grouping_ids)}")
    logger.info(f"  - Missing grouping (ERROR): {len(missing_grouping)}")
    logger.info(f"  - Zero viral hits (OK): {len(zero_vv)}")
    if missing_grouping:
        logger.error(f"Data integrity error: {len(missing_grouping)} samples have virus hits but missing grouping")
        for sample_id in sorted(missing_grouping):
            logger.error(f"  - Sample {sample_id} has virus hits but no grouping information")
    if zero_vv:
        logger.info(f"Found {len(zero_vv)} samples with zero viral hits (will be excluded from downstream analysis)")
        for sample_id in sorted(list(zero_vv)[:5]):  # Log first 5
            logger.info(f"  - Sample {sample_id} has no virus hits (zero viral hits)")
        if len(zero_vv) > 5:
            logger.info(f"  - ... and {len(zero_vv) - 5} more")
    return missing_grouping, zero_vv

#=======================================================================
# File processing functions
#=======================================================================

def filter_grouping_file(grouping_file: str, join_field: str, zero_vv_ids: set[str], validated_output: str) -> None:
    """
    Create filtered grouping file containing only samples with viral hits.
    
    This function reads the original grouping file and creates a new file
    that excludes samples with zero viral hits, preserving all other data.
    
    Args:
        grouping_file: Path to original grouping file (may be gzipped)
        join_field: Field name to filter on (must exist in header)
        zero_vv_ids: Set of sample IDs to exclude (zero viral hits)
        validated_output: Path to write filtered grouping file (gzipped if .gz extension)
    
    Returns:
        None
    """
    input_handle = open_by_suffix(grouping_file, 'r')
    output_handle = open_by_suffix(validated_output, 'w')
    try:
        header_line = input_handle.readline().strip()
        header = header_line.split('\t')
        join_field_index = header.index(join_field)
        output_handle.write(header_line + '\n')
        rows_kept = 0
        rows_excluded = 0
        for line in input_handle:
            line = line.strip()
            fields = line.split('\t')
            sample_id = fields[join_field_index]
            if sample_id in zero_vv_ids:
                rows_excluded += 1
                logger.debug(f"Excluding sample {sample_id} (zero viral hits)")
            else:
                output_handle.write(line + '\n')
                rows_kept += 1
        logger.info(f"Filtered grouping file: {rows_kept} rows kept, {rows_excluded} rows excluded")
    finally:
        input_handle.close()
        output_handle.close()

def write_zero_vv_log(zero_vv_ids: set[str], samples_without_vv_hits: str) -> None:
    """Write zero viral hits log file.
    
    Creates a TSV log file documenting samples that have no viral hits.
    
    Args:
        zero_vv_ids: Set of sample IDs with zero viral hits
        samples_without_vv_hits: Path to write the log file
    
    Returns:
        None
    """
    with open(samples_without_vv_hits, 'w') as f:
        f.write("sample_id\n")
        for sample_id in sorted(zero_vv_ids):
            f.write(f"{sample_id}\n")
    logger.info(f"Wrote zero viral hits log: {len(zero_vv_ids)} samples to {samples_without_vv_hits}")

#=======================================================================
# Main function
#=======================================================================

def main() -> None:
    """Main function implementing directional validation logic.
    
    This function orchestrates the validation process:
    1. Reads sample IDs from both virus hits and grouping files
    2. Performs directional validation to check data integrity
    3. Fails if samples with virus hits are missing grouping
    4. Filters grouping file to exclude zero-VV samples
    5. Creates log file documenting zero-VV samples
    
    Args:
        None (uses command-line arguments)
    
    Returns:
        None
    
    Raises:
        ValueError: If samples with virus hits are missing grouping information
    """
    args = parse_args()
    logger.info("Starting directional validation")
    logger.info(f"Virus hits file: {args.virus_hits_file}")
    logger.info(f"Grouping file: {args.grouping_file}")
    logger.info(f"Join field: {args.join_field}")
    # Step 1: Read sample IDs from both files
    virus_hits_ids, _ = read_sample_ids_from_file(args.virus_hits_file, args.join_field)
    grouping_ids, _ = read_sample_ids_from_file(args.grouping_file, args.join_field)
    # Step 2: Perform directional validation
    missing_grouping_ids, zero_vv_ids = validate_directional_integrity(virus_hits_ids, grouping_ids)
    # Step 3: Check for data integrity errors (FAIL if samples with VVs missing grouping)
    if missing_grouping_ids:
        error_msg = f"Validation failed: {len(missing_grouping_ids)} samples with virus hits are missing grouping information"
        logger.error(error_msg)
        raise ValueError(error_msg)
    # Step 4: Filter grouping file to remove zero-VV samples
    filter_grouping_file(args.grouping_file, args.join_field, zero_vv_ids, args.validated_output)
    # Step 5: Write zero-VV log for audit trail
    write_zero_vv_log(zero_vv_ids, args.samples_without_vv_hits)
    logger.info("Directional validation completed successfully")
    logger.info(f"Validated grouping file: {args.validated_output}")
    logger.info(f"Zero VV log file: {args.samples_without_vv_hits}")

if __name__ == "__main__":
    main()
