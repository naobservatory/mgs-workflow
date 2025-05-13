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
from dataclasses import dataclass
from collections import defaultdict

# Configure logging
class UTCFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime('%Y-%m-%d %H:%M:%S UTC')
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter('[%(asctime)s] %(message)s')
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)

# Define constants
TAXID_ROOT = 1

# Define Group class
@dataclass
class Group:
    group_id: str
    taxids: set[int]
    n_entries: int
    min_score: float
    max_score: float
    mean_score: float
    lca: int | None

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
    parent_to_children: dict[int, set[int]] = defaultdict(set)
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
                logger.warning(msg)
                path.append(TAXID_ROOT)
                break
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
    logger.debug(f"Finding pairwise LCA for taxids: {taxid1} and {taxid2}")
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
    logger.debug(f"Found LCA for taxids {taxid1} and {taxid2}: {lca}")
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
    logger.debug(f"Finding LCA for taxids: {taxids}")
    taxids_start = taxids.copy()
    assert isinstance(taxids, set), f"Taxids must be a set, got {type(taxids)} ({taxids})."
    assert len(taxids) > 0, "Taxids must be non-empty."
    # If only one taxid, return it as LCA
    if len(taxids) == 1:
        taxid = taxids.pop()
        logger.debug(f"Only one taxid, returning it as LCA: {taxid}")
        return taxid, path_cache
    # If any taxid is not in the dictionary, raise a warning and return the root
    for taxid in taxids:
        if taxid not in child_to_parent:
            logger.warning(f"Taxid {taxid} not found in child_to_parent dictionary; returning root.")
            return TAXID_ROOT, path_cache
    # Otherwise, find iterated pairwise LCA
    lca, path_cache = find_lca_pair(taxids.pop(), taxids.pop(), child_to_parent, path_cache)
    while taxids:
        lca, path_cache = find_lca_pair(lca, taxids.pop(), child_to_parent, path_cache)
    # Return the LCA and the updated path cache
    logger.debug(f"Found LCA for taxids {taxids_start}: {lca}")
    return lca, path_cache

def new_group(
        group_id: str,
        taxid: int,
        score: float,
        ) -> Group:
    """
    Create a new Group object from minimal input information.
    Args:
        group_id (str): Group ID.
        taxid (int): Taxid.
        score (float): Score.
    Returns:
        Group: New Group object.
    """
    return Group(
        group_id=group_id,
        taxids=set([taxid]),
        n_entries=1,
        min_score=score,
        max_score=score,
        mean_score=score,
        lca=None,
    )

def process_input_line(
        fields: list[str],
        group_idx: int,
        taxid_idx: int,
        score_idx: int,
        group_info: Group | None,
        ) -> Group:
    """
    Process a single line of an input TSV and return an updated Group object
    containing information about the corresponding entry group.
    Args:
        fields (list[str]): List of fields from an input line.
        group_idx (int): Index of the group field in the input line.
        taxid_idx (int): Index of the taxid field in the input line.
        score_idx (int): Index of the score field in the input line.
        group_info (Group | None): Group object containing information about the
            corresponding entry group.
    Returns:
        Group: Updated Group object containing information about the
            corresponding entry group.
    """
    logger.debug(f"Processing input line: {fields}")
    logger.debug(f"Group info: {group_info}")
    # Get group ID and confirm match with Group object
    group_id = fields[group_idx]
    if group_info is not None:
        assert group_id == group_info.group_id, \
            f"Group ID mismatch: {group_id} != {group_info.group_id}"
    # Get taxid and score
    taxid = int(fields[taxid_idx])
    score = float(fields[score_idx])
    # If no Group object, create one
    if group_info is None:
        return new_group(group_id, taxid, score)
    # Otherwise, update Group object
    min_score = min(group_info.min_score, score)
    max_score = max(group_info.max_score, score)
    score_sum = group_info.mean_score * group_info.n_entries + score
    mean_score = score_sum / (group_info.n_entries + 1)
    # Create a new set with all the existing taxids plus the new one
    updated_taxids = group_info.taxids.copy()
    updated_taxids.add(taxid)
    return Group(
        group_id=group_id,
        taxids=updated_taxids,
        n_entries=group_info.n_entries + 1,
        min_score=min_score,
        max_score=max_score,
        mean_score=mean_score,
        lca=group_info.lca,
    )

def process_group_to_output(
        group_info: Group,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
    ) -> tuple[str, dict[int, list[int]]]:
    """
    Process a Group object and return a string for output.
    Args:
        group_info (Group): Group object to process.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
    Returns:
        tuple[str, dict[int, list[int]]]: Tuple containing the output line
            and the updated path cache.
    """
    logger.debug(f"Processing group: {group_info}")
    # Get LCA of taxids in group
    lca, path_cache = find_lca_set(group_info.taxids, child_to_parent, path_cache)
    group_info.lca = lca
    # Prepare output line
    out_fields = [
        group_info.group_id,
        str(lca),
        str(group_info.n_entries),
        str(group_info.min_score),
        str(group_info.max_score),
        str(group_info.mean_score),
    ]
    out_line = "\t".join(out_fields) + "\n"
    # Return output line and updated path cache
    return out_line, path_cache

def get_output_header() -> str:
    """
    Get the header for the output TSV.
    Returns:
        str: Header for the output TSV.
    """
    fields = ["group_id", "lca", "n_entries", "min_score", "max_score", "mean_score"]
    return "\t".join(fields) + "\n"

def parse_input_header(
        header_fields: list[str],
        group_field: str,
        taxid_field: str,
        score_field: str,
        ) -> tuple[int, int, int]:
    """
    Parse the header of an input TSV and return the indices of the group,
    taxid, and score fields.
    Args:
        header_fields (list[str]): List of fields from the header of an input TSV.
        group_field (str): Column header for group field.
        taxid_field (str): Column header for taxid field.
        score_field (str): Column header for score field.
    Returns:
        tuple[int, int, int]: Tuple containing the indices of the group,
            taxid, and score fields.
    """
    # Confirm fields are present
    assert group_field in header_fields, f"Group field {group_field} not found in header."
    assert taxid_field in header_fields, f"Taxid field {taxid_field} not found in header."
    assert score_field in header_fields, f"Score field {score_field} not found in header."
    # Return indices
    group_idx = header_fields.index(group_field)
    taxid_idx = header_fields.index(taxid_field)
    score_idx = header_fields.index(score_field)
    return group_idx, taxid_idx, score_idx

def parse_input_tsv(
        input_path: str,
        output_path: str,
        group_field: str,
        taxid_field: str,
        score_field: str,
        child_to_parent: dict[int, int],
        ) -> None:
    """
    Iterate linewise over input TSV, calculating the LCA for each group
    of entries and writing the result to the output file.
    Args:
        input_path (str): Path to input TSV.
        output_path (str): Path to output TSV.
        group_field (str): Column header for group field.
        taxid_field (str): Column header for taxid field.
        score_field (str): Column header for score field.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
    """
    with open_by_suffix(input_path) as inf, open_by_suffix(output_path, "w") as outf:
        # Write header to output file
        outf.write(get_output_header())
        # Read header from input file
        header = inf.readline().strip().split("\t")
        # Get indices of fields
        group_idx, taxid_idx, score_idx = parse_input_header(
            header, group_field, taxid_field, score_field)
        def parse_line(fields: list[str], group_info: Group | None) -> Group:
            return process_input_line(fields, group_idx, taxid_idx, score_idx, group_info)
        # Process first line
        fields = inf.readline().strip().split("\t")
        group_id = fields[group_idx]
        group_info = parse_line(fields, None)
        # Iterate over input file
        n_entries = 1
        n_groups = 1
        path_cache = {}
        for line in inf:
            # Parse fields
            fields = line.strip().split("\t")
            # Get group ID and check sorting
            group_id_new = fields[group_idx]
            assert group_id_new >= group_id, \
                f"Group ID out of order: {group_id_new} < {group_id}"
            # If group ID is unchanged, process line into existing Group object
            if group_id_new == group_id:
                group_info = parse_line(fields, group_info)
            # Otherwise, write previous Group object to output and start new one
            else:
                out_line, path_cache = process_group_to_output(
                    group_info, child_to_parent, path_cache)
                outf.write(out_line)
                group_info = parse_line(fields, None)
                group_id = group_id_new
                n_groups += 1
            n_entries += 1
        # Write last Group object to output
        out_line, path_cache = process_group_to_output(
            group_info, child_to_parent, path_cache)
        outf.write(out_line)
        logger.info(f"Processed {n_entries} entries in {n_groups} groups.")

#=======================================================================
# Main function
#=======================================================================

def main() -> None:
    logger.info("Initializing script.")
    # Parse arguments
    logger.info("Parsing arguments.")
    args = parse_args()
    # Import taxonomy DB and process into dictionaries
    logger.info("Parsing taxonomy DB.")
    child_to_parent, parent_to_children = parse_taxonomy_db(args.db)
    logger.info(f"Parsed taxonomy information for {len(child_to_parent)} taxids.")
    # Parse input TSV and write LCA information to output TSV
    logger.info("Parsing input TSV.")
    parse_input_tsv(args.input, args.output, args.group, args.taxid, args.score, child_to_parent)
    # Log completion
    logger.info("Script complete.")

if __name__ == "__main__":
    main()