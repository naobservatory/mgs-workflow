#!/usr/bin/env python

"""
Given a sorted TSV containing group, taxid, and score columns,
return a TSV with a single row per group containing the lowest common ancestor
taxid across all taxids in the group.
"""

#=======================================================================
# Preamble
#=======================================================================

# Import libraries
import logging
import argparse
from datetime import datetime, timezone
import gzip
import bz2

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

# Define constants
TAXID_ROOT = 1

#=======================================================================
# Auxiliary functions
#=======================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    # Create parser
    desc = "Given a sorted TSV containing group, taxid, and score columns, " \
           "return a TSV with a single row per group containing the lowest " \
           "common ancestor taxid across all taxids in the group."
    parser = argparse.ArgumentParser(description=desc)
    # Add arguments
    parser.add_argument("--input", "-i", help="Path to input TSV.")
    parser.add_argument("--output", "-o", help="Path to output TSV.")
    parser.add_argument("--db", "-d", help="Path to taxonomy DB (raw NCBI nodes.dmp file).")
    parser.add_argument("--group", "-g", help="Column header for group field.")
    parser.add_argument("--taxid", "-t", help="Column header for taxid field.")
    parser.add_argument("--score", "-s", help="Column header for score field.")
    # Return parsed arguments
    return parser.parse_args()

def open_by_suffix(filename, mode="r", debug=False):
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

def parse_taxonomy_db(path: str) -> dict:
    """Parse taxonomy DB into two dictionaries: one mapping each taxid
    to its parent taxid, and one mapping each taxid to its children taxids."""
    # Define dictionaries
    child_to_parent: dict[int, int] = {}
    parent_to_children: dict[int, set[int]] = {}
    # Read file line by line and parse into dictionaries
    with open_by_suffix(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            taxid = int(fields[0])
            parent_taxid = int(fields[2])
            child_to_parent[taxid] = parent_taxid
            parent_to_children[parent_taxid].add(taxid)
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
    # Check input type
    assert isinstance(taxid, int), "Taxid must be an integer."
    # If path is already cached, return it
    if taxid in path_cache:
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
                logger.warning(msg)
                path.append(TAXID_ROOT)
                break
            path.append(parent)
    # Check path includes taxid and root
    assert path[0] == taxid, "Path does not start with taxid."
    assert path[-1] == TAXID_ROOT, "Path does not end with root."
    # Add path to cache
    path_cache[taxid] = path
    # Return the path and cache
    return path, path_cache

def find_lca_paths(
        path1: list[int],
        path2: list[int],
        ) -> int:
    """
    Given two paths from a taxid to the root, find the lowest common ancestor
    shared by both paths.
    Args:
        path1 (list[int]): Path from first taxid to root.
        path2 (list[int]): Path from second taxid to root.
    Returns:
        int: The LCA taxid.
    """
    # If paths are identical, return the first taxid
    if path1 == path2:
        return path1[0]
    # Otherwise, walk down paths from root until they diverge
    lca = TAXID_ROOT
    while path1 and path2:
        ancestor1 = path1.pop()
        ancestor2 = path2.pop()
        if ancestor1 != ancestor2:
            break
        lca = ancestor1
    return lca

def find_lca_pair(
        taxid1: int,
        taxid2: int,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        ) -> tuple[int, dict[int, list[int]]]:
    """
    Find the lowest common ancestor (LCA) of two taxids.
    Args:
        taxid1 (int): First taxid.
        taxid2 (int): Second taxid.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
    Returns:
        tuple[int, dict[int, list[int]]]: Tuple containing the LCA taxid
            and the updated path cache.
    """
    # Check types
    assert isinstance(taxid1, int), "Taxid 1 must be an integer."
    assert isinstance(taxid2, int), "Taxid 2 must be an integer."
    # If taxids are the same, return that value as LCA
    if taxid1 == taxid2:
        return taxid1, path_cache
    # Get paths to root
    path1, path_cache = path_to_root(taxid1, child_to_parent, path_cache)
    path2, path_cache = path_to_root(taxid2, child_to_parent, path_cache)
    # Find LCA by walking down paths until they diverge
    lca = find_lca_paths(path1, path2)
    return lca, path_cache

def find_lca_set(
        taxids: set[int],
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        ) -> tuple[int, dict[int, list[int]]]:
    """
    Find the lowest common ancestor (LCA) of a set of taxids.
    Args:
        taxids (set[int]): Set of taxids to find the LCA for.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
    Returns:
        tuple[int, dict[int, list[int]]]: Tuple containing the LCA taxid
            and the updated path cache.
    """
    assert isinstance(taxids, set), "Taxids must be a set."
    assert len(taxids) > 0, "Taxids must be non-empty."
    # If only one taxid, return it as LCA
    if len(taxids) == 1:
        return taxids.pop(), path_cache
    # If any taxid is not in the dictionary, raise a warning and return the root
    for taxid in taxids:
        if taxid not in child_to_parent:
            logger.warning(f"Taxid {taxid} not found in child_to_parent dictionary.")
            return TAXID_ROOT, path_cache
    # Otherwise, find iterated pairwise LCA
    lca, path_cache = find_lca_pair(taxids.pop(), taxids.pop(), child_to_parent, path_cache)
    while taxids:
        lca, path_cache = find_lca_pair(lca, taxids.pop(), child_to_parent, path_cache)
    # Return the LCA and the updated path cache
    return lca, path_cache

#=======================================================================
# Main function
#=======================================================================

def main() -> None:
    logger.info("Initializing script.")
    # Parse arguments
    args = parse_args()
    # Import taxonomy DB and process into dictionaries
    child_to_parent, parent_to_children = parse_taxonomy_db(args.db)

    logger.info("Script complete.")

if __name__ == "__main__":
    main()