#!/usr/bin/env python

"""
Write a new grouping file with only samples observed in the hits file, and write a list of samples without viral hits.

Usage:
    validate_grouping.py virus_hits.tsv grouping.tsv validated_grouping.tsv.gz samples_without_vv_hits.tsv
"""

#=======================================================================
# Preamble
#=======================================================================

# Import libraries
import argparse
import logging
import gzip
import bz2
from datetime import datetime, timezone
from typing import TextIO
import os

# Configure logging
class UTCFormatter(logging.Formatter):
    """Custom formatter that outputs timestamps in UTC format.
    
    This formatter ensures all log timestamps are in UTC for consistency
    across different timezones and environments.
    """
    
    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        """Format the time in UTC.
        
        Args:
            record (logging.LogRecord): LogRecord instance containing the event data
            datefmt (str | None): Date format string (currently unused, UTC format is hardcoded)
        
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
        filename (str): Path to the file to open
        mode (str): File open mode (default: "r")
    
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
    
    Returns:
        argparse.Namespace: Parsed command-line arguments containing:
            - virus_hits_file: Path to virus hits TSV file
            - grouping_file: Path to grouping TSV file
            - validated_output: Path to validated grouping output file
            - samples_without_vv_hits: Path to samples without viral hits log file
    """
    parser = argparse.ArgumentParser(description="Validate grouping file against virus hits file")
    parser.add_argument("virus_hits_file", help="Path to virus hits TSV file (must contain a 'sample' column)")
    parser.add_argument("grouping_file", help="Path to grouping TSV file (columns: group	sample)")
    parser.add_argument("validated_output", help="Path to validated grouping output file (.tsv or .tsv.gz)")
    parser.add_argument("samples_without_vv_hits", help="Path to samples-without-viral-hits output (.tsv or .tsv.gz)")
    return parser.parse_args()

#=======================================================================
# Core logic (no per-sample counts needed)
#=======================================================================

def _read_grouping(grouping_path: str) -> dict[str, str]:
    """Load the entire grouping TSV, finding 'group' and 'sample' columns by header.
    
    Args:
        grouping_path (str): Path to the grouping TSV file
    
    Returns:
        dict[str, str]: Mapping from sample ID to group ID
    
    Raises:
        ValueError: If 'group' or 'sample' columns not found, duplicate samples found, or file is empty
    """
    sample_to_group: dict[str, str] = {}
    header_idx = None
    
    with open_by_suffix(grouping_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header_idx is None:
                header_idx = {name.lower(): i for i, name in enumerate(parts)}
                if 'group' not in header_idx or 'sample' not in header_idx:
                    raise ValueError("Grouping file header must contain 'group' and 'sample' columns")
                continue
            group = parts[header_idx['group']]
            sample = parts[header_idx['sample']]
            if sample in sample_to_group:
                raise ValueError(f"Duplicate sample in grouping: {sample}")
            sample_to_group[sample] = group
    if not sample_to_group:
        raise ValueError("Grouping file is empty (no data rows)")
    logger.info(f"Loaded grouping with {len(sample_to_group)} samples")
    return sample_to_group


def _stream_hits(virus_hits_path: str, sample_to_group: dict[str, str]) -> set[str]:
    """Stream virus_hits TSV, return set of sample IDs that appeared in virus_hits.
    
    Args:
        virus_hits_path (str): Path to the virus hits TSV file
        sample_to_group (dict[str, str]): Mapping from sample ID to group ID
    
    Returns:
        set[str]: Set of sample IDs that appeared in virus_hits
    """
    seen: set[str] = set()
    header_idx = None
    total_rows = 0
    with open_by_suffix(virus_hits_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header_idx is None:
                header_idx = {name: i for i, name in enumerate(parts)}
                if 'sample' not in header_idx:
                    raise ValueError("Virus hits header must contain a 'sample' column")
                continue
            total_rows += 1
            sample = parts[header_idx['sample']]
            seen.add(sample)
    logger.info(f"Streamed {total_rows} virus-hit rows; unique samples observed: {len(seen)}")
    missing = seen - set(sample_to_group.keys())
    if missing:
        sorted_missing = sorted(missing)
        logger.debug(
            f"Found {len(missing)} samples present in virus_hits but MISSING from grouping: "
            f"{', '.join(sorted_missing[:10])}{'...' if len(missing) > 10 else ''}"
        )
    return seen

def _write_outputs(
    sample_to_group: dict[str, str],
    seen_in_hits: set[str],
    out_new_grouping_path: str,
    out_no_vv_samples_path: str,
) -> None:
    """Write samples with no viral hits (grouping - hits) and new grouping (intersection).
    
    Args:
        sample_to_group (dict[str, str]): Mapping from sample ID to group ID
        seen_in_hits (set[str]): Set of sample IDs observed in virus hits
        out_new_grouping_path (str): Path to write validated grouping file
        out_no_vv_samples_path (str): Path to write samples with no viral hits
    """
    grouping_samples = set(sample_to_group.keys())
    samples_no_viral_hits = sorted(grouping_samples - seen_in_hits)
    with open_by_suffix(out_no_vv_samples_path, 'w') as out_unused:
        for s in samples_no_viral_hits:
            out_unused.write(f"{s}\n")
    logger.info(f"Wrote {len(samples_no_viral_hits)} samples with no viral hits to {out_no_vv_samples_path}")
    samples_with_viral_hits = sorted(grouping_samples & seen_in_hits)
    with open_by_suffix(out_new_grouping_path, 'w') as out_group:
        out_group.write("group\tsample\n")
        for s in samples_with_viral_hits:
            out_group.write(f"{sample_to_group[s]}\t{s}\n")
    logger.info(f"Wrote validated grouping with {len(samples_with_viral_hits)} samples to {out_new_grouping_path}")


def validate_grouping(
    grouping_path: str,
    virus_hits_path: str,
    out_no_vv_samples_path: str,
    out_new_grouping_path: str,
) -> None:
    """Update grouping file to only include samples observed in the viral hits file.
    
    Args:
        grouping_path (str): Path to the grouping TSV file
        virus_hits_path (str): Path to the virus hits TSV file
        out_no_vv_samples_path (str): Path to write samples with no viral hits
        out_new_grouping_path (str): Path to write validated grouping file
    
    Raises:
        ValueError: If virus hits contain samples missing from grouping
    """
    # 1) Load grouping
    sample_to_group = _read_grouping(grouping_path)
    # 2) Stream hits and validate that samples are subset of grouping
    seen_in_hits = _stream_hits(virus_hits_path, sample_to_group)
    # 3) Write samples with no viral hits to out_no_vv_samples_path and new grouping to out_new_grouping_path
    _write_outputs(sample_to_group, seen_in_hits, out_new_grouping_path, out_no_vv_samples_path)


#=======================================================================
# Entrypoint
#=======================================================================

def main() -> None:
    """Parse command-line arguments and call validate_grouping().

    Raises:
        ValueError: If validation fails due to missing samples in grouping
    """
    args = parse_args()
    logger.info("Starting validation")
    try:
        validate_grouping(
            grouping_path=args.grouping_file,
            virus_hits_path=args.virus_hits_file,
            out_no_vv_samples_path=args.samples_without_vv_hits,
            out_new_grouping_path=args.validated_output,
        )
    except ValueError as e:
        logger.error(str(e))
        raise
    logger.info("Validation complete")


if __name__ == "__main__":
    main()

