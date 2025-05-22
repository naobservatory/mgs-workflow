#!/usr/bin/env python

"""
Given a TSV with two taxid columns, compute the vertical taxonomic distance
(number of parent-child steps) between the two taxids for each row. Rows
for which the two taxids are the same are given a distance of 0; rows for which
one taxid is not an ancestor of the other are given a distance of NA.
"""

#=======================================================================
# Import libraries
#=======================================================================

import logging
import argparse
from datetime import datetime, timezone
import time
import gzip
import bz2
from dataclasses import dataclass
from collections import defaultdict
from typing import TextIO

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
# Define constants
#=======================================================================

TAXID_ROOT = 1

#=======================================================================
# I/O functions
#=======================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    # Create parser
    desc = (
        "Given a TSV with two taxid columns, compute the vertical taxonomic "
        "distance (number of parent-child steps) between the two taxids for "
        "each row. Rows for which the two taxids are the same are given a "
        "distance of 0; rows for which one taxid is not an ancestor of the "
        "other are given a distance of NA."
    )
    parser = argparse.ArgumentParser(description=desc)
    # Add arguments
    parser.add_argument("--input", "-i", help="Path to input TSV.")
    parser.add_argument("--output", "-o", help="Path to output TSV.")
    parser.add_argument("--taxid-field-1", "-t1", help="Column header for first taxid field.")
    parser.add_argument("--taxid-field-2", "-t2", help="Column header for second taxid field.")
    parser.add_argument("--distance-field", "-d", help="Column header for new distance field.")
    parser.add_argument("--nodes-db", "-n", help="Path to taxonomy nodes DB (raw NCBI nodes.dmp file).")
    # Return parsed arguments
    return parser.parse_args()

def open_by_suffix(filename, mode="r", debug=False):
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
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

def parse_header(header_line: str, fields: list[str]) -> tuple[list[str], dict[str, int]]:
    """
    Parse a TSV header line into a list of fields, and compute the
    indices of a set of expected fields.
    Args:
        header_line (str): Header line of a TSV.
        fields (list[str]): List of field names to check for in the header line.
    Returns:
        tuple[list[str], dict[str, int]]: Tuple containing the list of fields
            from the header line and a dictionary mapping field names to their
            indices.
    """
    # Check for empty file
    header_line_stripped = header_line.strip()
    if not header_line_stripped:
        raise ValueError("Header line is empty: no fields to parse.")
    # Split header line into fields
    header_fields = header_line_stripped.split("\t")
    # Check that all expected fields are present and get their indices
    indices = {}
    for field in fields:
        if field not in header_fields:
            raise ValueError(f"Field not found in header: {field}")
        indices[field] = header_fields.index(field)
    # Return fields and indices
    return header_fields, indices

#=======================================================================
# Taxonomy functions
#=======================================================================

def parse_taxid(taxid_str: str) -> int|None:
    """Parse a taxid string into an integer."""
    try:
        return int(taxid_str)
    except ValueError:
        return None

def parse_nodes_db(path: str
                      ) -> tuple[dict[int, int], dict[int, set[int]]]:
    """
    Parse taxonomy DB into two dictionaries: one mapping each taxid
    to its parent taxid, and one mapping each taxid to its children taxids.
    Args:
        path (str): Path to taxonomy DB.
    Returns:
        tuple[dict[int, int], dict[int, set[int]]]: Tuple containing the
            child-to-parent dictionary and the parent-to-children dictionary.
    """
    # Define dictionaries
    child_to_parent: dict[int, int] = {}
    parent_to_children: dict[int, set[int]] = defaultdict(set)
    # Read file line by line and parse into dictionaries
    with open_by_suffix(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            # Parse taxids strictly (not tolerating non-integer strings)
            taxid = int(fields[0])
            parent_taxid = int(fields[2])
            child_to_parent[taxid] = parent_taxid
            parent_to_children[parent_taxid].add(taxid)
    # Check that DB contains root as the topmost taxid
    assert TAXID_ROOT in child_to_parent and TAXID_ROOT in parent_to_children, \
        "Taxonomy DB does not contain root."
    assert child_to_parent[TAXID_ROOT] == TAXID_ROOT, "Root taxid has a parent." # NCBI file has root as a child of itself
    assert TAXID_ROOT in parent_to_children[TAXID_ROOT], "Root taxid must be its own child."
    # Return dictionaries
    return child_to_parent, parent_to_children

def path_to_root(
        taxid: int,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        ) -> tuple[list[int], dict[int, list[int]]]:
    """
    Find the path from a taxid to the root of the taxonomy tree.
    Args:
        taxid (int): The starting taxid.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
    Returns:
        tuple[list[int], dict[int, list[int]]]: Tuple containing the path to the root
            and the updated path cache.
    """
    logger.debug(f"Finding path to root for taxid: {taxid}")
    # Check input type
    assert isinstance(taxid, int), "Taxid must be an integer."
    # If path is already cached, return it
    if taxid in path_cache:
        logger.debug(f"Path to root for target taxid {taxid} already cached: {path_cache[taxid]}")
        return path_cache[taxid], path_cache
    # Initialize path
    path: list[int] = [taxid]
    # If taxid is not in child_to_parent, raise a warning and return the root
    if taxid not in child_to_parent:
        logger.warning(f"Taxid {taxid} not found in child_to_parent dictionary.")
        path.append(TAXID_ROOT)
    # Otherwise, traverse up the tree until the root is reached
    else:
        while path[-1] != TAXID_ROOT:
            parent = child_to_parent[path[-1]]
            # Should not encounter loops before reaching the root
            if parent == path[-1]:
                msg = f"Taxid {taxid} has a self-loop in the child_to_parent dictionary."
                logger.error(msg)
                raise ValueError(msg)
            # Check if path for parent is already cached
            if parent in path_cache:
                logger.debug(f"Path to root for ancestor taxid {parent} already cached: {path_cache[parent]}")
                path.extend(path_cache[parent])
            # Otherwise, add parent to path
            else:
                path.append(parent)
    # Check path includes taxid and root
    assert path[0] == taxid, "Path does not start with taxid."
    assert path[-1] == TAXID_ROOT, "Path does not end with root."
    # Add paths to cache for target taxid and all its parents
    for i in range(len(path)):
        # Each path is a suffix of the previous path
        path_cache[path[i]] = path[i:]
    # Return the path and cache
    logger.debug(f"Path to root for target taxid {taxid}: {path}")
    return path, path_cache

def compute_taxonomic_distance(
        taxid_1: int,
        taxid_2: int,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        ) -> tuple[int|None, dict[int, list[int]]]:
    """
    Compute the taxonomic distance between two taxids as the number of
    parent-child steps between them. Distance is negative if taxid_1 is an
    ancestor of taxid_2 and positive if taxid_2 is an ancestor of taxid_1.
    Args:
        taxid_1 (int): The first taxid.
        taxid_2 (int): The second taxid.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
    Returns:
        tuple[int, dict[int, list[int]]]: Tuple containing the taxonomic distance
            and the updated path cache.
    """
    # If taxids are the same, return 0
    if taxid_1 == taxid_2:
        return 0, path_cache
    if taxid_1 is None or taxid_2 is None:
        return None, path_cache
    # Get paths to root for both taxids (starting from taxid itself)
    path_1, path_cache = path_to_root(taxid_1, child_to_parent, path_cache)
    logger.debug(f"Path to root for taxid {taxid_1}: {path_1}")
    path_2, path_cache = path_to_root(taxid_2, child_to_parent, path_cache)
    logger.debug(f"Path to root for taxid {taxid_2}: {path_2}")
    # Compute taxonomic distance
    # Distance = d1 - d2, where dN is the number of steps from taxid N to the root
    # Therefore, if taxid_1 is an ancestor of taxid_2, the distance is negative
    if taxid_1 in path_2:
        distance = -path_2.index(taxid_1)
    elif taxid_2 in path_1:
        distance = path_1.index(taxid_2)
    else:
        distance = None
    logger.debug(f"Taxonomic distance between taxids {taxid_1} and {taxid_2}: {distance}")
    # Return distance and path cache
    return distance, path_cache

#=======================================================================
# Functions for processing input and output
#=======================================================================

def process_input_to_output(
        input_path: str,
        output_path: str,
        field_names: dict[str, str],
        child_to_parent: dict[int, int],
        ) -> None:
    """
    Iterate linewise over input TSV, computing the taxonomic distance
    between the two taxids for each group of entries and writing the result
    to the output file.
    Args:
        input_path (str): Path to input TSV.
        output_path (str): Path to output TSV.
        field_names (dict[str, str]): Dictionary containing taxid and distance field names.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
    """
    with open_by_suffix(input_path) as inf, open_by_suffix(output_path, "w") as outf:
        # Read and handle input header
        logger.info("Parsing input header.")
        fields_to_check = [field_names["taxid_1"], field_names["taxid_2"]]
        header_line = inf.readline().strip()
        header_fields, indices = parse_header(header_line, fields_to_check)
        logger.info(f"Parsed input header: {header_fields}")
        if field_names["distance"] in header_fields:
            msg = (
                "Distance field already present in input header: "
                f"{field_names['distance']}. "
            )
            raise ValueError(msg)
        indices[field_names["distance"]] = len(header_fields)
        logger.info(f"Indices of target fields: {indices}")
        # Write output header
        header_fields_out = header_fields + [field_names["distance"]]
        outf.write(join_line(header_fields_out))
        # Process rest of input file
        path_cache = {}
        n_entries = 0
        for line in inf:
            # If line is empty, break
            if not line.strip():
                break 
            # Parse line into fields and extract taxids (parsing non-integer strings as None)
            fields = line.strip().split("\t")
            taxid_1 = parse_taxid(fields[indices[field_names["taxid_1"]]])
            taxid_2 = parse_taxid(fields[indices[field_names["taxid_2"]]])
            # Compute taxonomic distance
            distance, path_cache = compute_taxonomic_distance(taxid_1, taxid_2,
                                                              child_to_parent,
                                                              path_cache)
            # Write output line
            distance_out = [str(distance)] if distance is not None else ["NA"]
            fields_out = fields + distance_out
            outf.write(join_line(fields_out))
            n_entries += 1
        logger.info(f"Processed {n_entries} entries.")

#=======================================================================
# Main function
#=======================================================================

def main() -> None:
    logger.info("Initializing script.")
    start_time = time.time()
    # Parse arguments
    logger.info("Parsing arguments.")
    args = parse_args()
    logger.info(f"Arguments: {args}")
    # Import taxonomy DB and process into dictionaries
    logger.info("Parsing taxonomy DB.")
    child_to_parent, parent_to_children = parse_nodes_db(args.nodes_db)
    logger.info(f"Parsed taxonomy information for {len(child_to_parent)} taxids.")
    logger.debug(f"Child-to-parent dictionary: {child_to_parent}")
    # Prepare fields dict
    fields = {"taxid_1": args.taxid_field_1,
              "taxid_2": args.taxid_field_2,
              "distance": args.distance_field}
    # Parse input TSV and compute taxonomic distances
    logger.info("Parsing input TSV and computing taxonomic distances.")
    process_input_to_output(args.input, args.output, fields, child_to_parent)
    # Log completion
    logger.info("Script completed successfully.")
    end_time = time.time()
    logger.info(f"Total time elapsed: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()