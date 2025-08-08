#!/usr/bin/env python

"""
Validate grouping file against virus hits file using directional validation.

This script implements asymmetric validation:
1. ENFORCE: Every sample in virus_hits MUST have corresponding grouping (data integrity)
2. ALLOW: Samples in grouping without virus hits (zero viral variants - filter out)

Usage:
    validate_grouping.py virus_hits.tsv grouping.tsv join_field validated_grouping.tsv.gz zero_vv_log.tsv
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

def open_by_suffix(filename, mode="r", debug=False):
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Validate grouping file against virus hits file")
    parser.add_argument("virus_hits_file", help="Path to virus hits TSV file")
    parser.add_argument("grouping_file", help="Path to grouping TSV file")
    parser.add_argument("join_field", help="Field to validate on (e.g., 'sample')")
    parser.add_argument("validated_output", help="Path to validated grouping output file (.tsv.gz)")
    parser.add_argument("zero_vv_log", help="Path to zero viral variants log file (.tsv)")
    return parser.parse_args()

#=======================================================================
# File reading functions
#=======================================================================

def read_sample_ids_from_file(file_path, join_field):
    """
    Read sample IDs from a TSV file.
    
    Returns:
        set: Set of sample IDs found in the file
        list: Header fields from the file
    """
    sample_ids = set()
    
    with open_by_suffix(file_path) as file_handle:
        # Read header
        header_line = file_handle.readline().strip()
        if not header_line:
            logger.warning(f"File {file_path} is empty")
            return set(), []
        
        header = header_line.split('\t')
        
        # Validate join field exists
        if join_field not in header:
            raise ValueError(f"Join field '{join_field}' not found in {file_path}. Available fields: {header}")
        
        join_field_index = header.index(join_field)
        logger.info(f"Reading sample IDs from {file_path}, join field '{join_field}' at index {join_field_index}")
        
        # Read data rows
        for line_num, line in enumerate(file_handle, start=2):
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) != len(header):
                logger.warning(f"Line {line_num} in {file_path} has {len(fields)} fields, expected {len(header)}")
                continue
            
            sample_id = fields[join_field_index]
            sample_ids.add(sample_id)
        
        logger.info(f"Found {len(sample_ids)} unique sample IDs in {file_path}")
        return sample_ids, header

#=======================================================================
# Validation functions
#=======================================================================

def validate_directional_integrity(virus_hits_ids, grouping_ids):
    """
    Perform directional validation.
    
    Args:
        virus_hits_ids: Set of sample IDs from virus hits file
        grouping_ids: Set of sample IDs from grouping file
    
    Returns:
        tuple: (missing_grouping_ids, zero_vv_ids)
            - missing_grouping_ids: samples in virus_hits but not in grouping (ERROR)
            - zero_vv_ids: samples in grouping but not in virus_hits (OK - zero VVs)
    """
    logger = logging.getLogger()
    
    # Find samples with virus hits but missing grouping (DATA INTEGRITY ERROR)
    missing_grouping = virus_hits_ids - grouping_ids
    
    # Find samples with grouping but no virus hits (ZERO VIRAL VARIANTS - OK)
    zero_vv = grouping_ids - virus_hits_ids
    
    logger.info(f"Validation results:")
    logger.info(f"  - Samples with virus hits: {len(virus_hits_ids)}")
    logger.info(f"  - Samples with grouping: {len(grouping_ids)}")
    logger.info(f"  - Missing grouping (ERROR): {len(missing_grouping)}")
    logger.info(f"  - Zero viral variants (OK): {len(zero_vv)}")
    
    if missing_grouping:
        logger.error(f"Data integrity error: {len(missing_grouping)} samples have virus hits but missing grouping")
        for sample_id in sorted(missing_grouping):
            logger.error(f"  - Sample {sample_id} has virus hits but no grouping information")
    
    if zero_vv:
        logger.info(f"Found {len(zero_vv)} samples with zero viral variants (will be excluded from downstream analysis)")
        for sample_id in sorted(list(zero_vv)[:5]):  # Log first 5
            logger.info(f"  - Sample {sample_id} has no virus hits (zero viral variants)")
        if len(zero_vv) > 5:
            logger.info(f"  - ... and {len(zero_vv) - 5} more")
    
    return missing_grouping, zero_vv

#=======================================================================
# File processing functions
#=======================================================================

def filter_grouping_file(grouping_file, join_field, zero_vv_ids, validated_output):
    """
    Create filtered grouping file containing only samples with viral variants.
    
    Args:
        grouping_file: Path to original grouping file
        join_field: Field name to filter on
        zero_vv_ids: Set of sample IDs to exclude (zero viral variants)
        validated_output: Path to write filtered grouping file
    """
    logger = logging.getLogger()
    
    # Handle compressed input
    if grouping_file.endswith('.gz'):
        input_handle = gzip.open(grouping_file, 'rt')
    else:
        input_handle = open(grouping_file, 'r')
    
    # Write compressed output if .gz extension, otherwise uncompressed
    if validated_output.endswith('.gz'):
        output_handle = gzip.open(validated_output, 'wt')
    else:
        output_handle = open(validated_output, 'w')
    
    try:
        # Read and write header
        header_line = input_handle.readline().strip()
        header = header_line.split('\t')
        join_field_index = header.index(join_field)
        
        output_handle.write(header_line + '\n')
        
        # Filter data rows
        rows_kept = 0
        rows_excluded = 0
        
        for line in input_handle:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            sample_id = fields[join_field_index]
            
            if sample_id in zero_vv_ids:
                rows_excluded += 1
                logger.debug(f"Excluding sample {sample_id} (zero viral variants)")
            else:
                output_handle.write(line + '\n')
                rows_kept += 1
        
        logger.info(f"Filtered grouping file: {rows_kept} rows kept, {rows_excluded} rows excluded")
        
    finally:
        input_handle.close()
        output_handle.close()

def write_zero_vv_log(zero_vv_ids, zero_vv_log):
    """Write zero viral variants log file."""
    logger = logging.getLogger()
    
    with open(zero_vv_log, 'w') as f:
        f.write("sample_id\treason\ttimestamp\n")
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        
        for sample_id in sorted(zero_vv_ids):
            f.write(f"{sample_id}\tzero_viral_variants\t{timestamp}\n")
    
    logger.info(f"Wrote zero viral variants log: {len(zero_vv_ids)} samples to {zero_vv_log}")

#=======================================================================
# Main function
#=======================================================================

def main():
    """Main function implementing directional validation logic."""
    args = parse_args()
    
    logger.info("Starting directional validation")
    logger.info(f"Virus hits file: {args.virus_hits_file}")
    logger.info(f"Grouping file: {args.grouping_file}")
    logger.info(f"Join field: {args.join_field}")
    
    # Step 1: Read sample IDs from both files
    virus_hits_ids, virus_hits_header = read_sample_ids_from_file(args.virus_hits_file, args.join_field)
    grouping_ids, grouping_header = read_sample_ids_from_file(args.grouping_file, args.join_field)
    
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
    write_zero_vv_log(zero_vv_ids, args.zero_vv_log)
    
    logger.info("Directional validation completed successfully")
    logger.info(f"Validated grouping file: {args.validated_output}")
    logger.info(f"Zero VV log file: {args.zero_vv_log}")

if __name__ == "__main__":
    main()
