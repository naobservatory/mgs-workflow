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
from typing import TextIO
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

#=======================================================================
# LCA functions
#=======================================================================

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
    path1_copy = path1.copy()
    path2_copy = path2.copy()
    lca = TAXID_ROOT
    while path1_copy and path2_copy:
        ancestor1 = path1_copy.pop()
        ancestor2 = path2_copy.pop()
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

#=======================================================================
# Subgroup class and functions
#=======================================================================

@dataclass
class Subgroup:
    """
    A subgroup of a Group object, tracking information about a subset of
    entries, including taxids, scores, and the top taxid and LCA.
    """
    taxids: set[int]
    n_entries: int
    min_score: float | None
    max_score: float | None
    mean_score: float | None
    top_taxid: int | None
    lca: int | None

def new_subgroup(
        taxid: int | None,
        score: float | None,
    ) -> Subgroup:
    """
    Create a new Subgroup object from minimal input information.
    Args:
        taxid (int | None): Taxid of entry, if any.
        score (float | None): Score of entry, if any.
    Returns:
        Subgroup: New Subgroup object.
    """
    return Subgroup(
        taxids = set([taxid]) if taxid is not None else set(),
        n_entries = 1 if taxid is not None else 0,
        min_score = score if score is not None else None,
        max_score = score if score is not None else None,
        mean_score = score if score is not None else None,
        top_taxid = taxid if taxid is not None else None,
        lca = taxid if taxid is not None else None,
    )

def update_subgroup(
        subgroup: Subgroup,
        taxid: int,
        score: float,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
    ) -> tuple[Subgroup, dict[int, list[int]]]:
    """
    Update a Subgroup object with new taxid and score information.
    Args:
        subgroup (Subgroup): Subgroup object to update.
        taxid (int): Taxid of new entry.
        score (float): Score of new entry.
    Returns:
        tuple[Subgroup, dict[int, list[int]]]: Tuple containing the updated
            Subgroup object and the updated path cache.
    """
    # Compute new LCA
    if subgroup.lca is None:
        lca, path_cache = taxid, path_cache
        min_score = score
        max_score = score
        mean_score = score
        top_taxid = taxid
    else:
        lca, path_cache = find_lca_pair(subgroup.lca, taxid, child_to_parent, path_cache)
        min_score = min(subgroup.min_score, score)
        max_score = max(subgroup.max_score, score)
        score_sum = subgroup.mean_score * subgroup.n_entries + score
        mean_score = score_sum / (subgroup.n_entries + 1)
        top_taxid = taxid if score > subgroup.max_score else subgroup.top_taxid
    # Update other attributes
    n_entries = subgroup.n_entries + 1
    taxids = subgroup.taxids.copy()
    taxids.add(taxid)
    # Return updated Subgroup object
    return Subgroup(
        taxids = taxids,
        n_entries = n_entries,
        min_score = min_score,
        max_score = max_score,
        mean_score = mean_score,
        top_taxid = top_taxid,
        lca = lca,
    ), path_cache

def output_subgroup(
        subgroup: Subgroup,
        ) -> list[str]:
    """
    Process a Subgroup object into a list of fields for output.
    """
    return [
        str(subgroup.lca),
        str(subgroup.n_entries),
        str(subgroup.top_taxid),
        str(subgroup.min_score),
        str(subgroup.max_score),
        str(subgroup.mean_score),
    ]

#=======================================================================
# Group class and functions
#=======================================================================

@dataclass
class Group:
    """
    A group of entries, including taxids, scores, and the top taxid and LCA,
    divided into natural and artificial subgroups.
    """
    group_id: str
    all: Subgroup
    natural: Subgroup
    artificial: Subgroup

def new_group(
        group_id: str,
        taxid: int,
        score: float,
        artificial_taxids: set[int],
        ) -> Group:
    """
    Create a new Group object from minimal input information.
    Args:
        group_id (str): Group ID.
        taxid (int): Taxid of first entry.
        score (float): Score of first entry.
        artificial_taxids (set[int]): Set of taxids that represent artificial
            or engineered sequences.
    Returns:
        Group: New Group object.
    """
    # Determine if first entry is artificial
    is_artificial = taxid in artificial_taxids
    # Initialize all, natural, and artificial subgroups
    subgroup_all = new_subgroup(taxid, score)
    if is_artificial:
        subgroup_natural = new_subgroup(None, None)
        subgroup_artificial = new_subgroup(taxid, score)
    else:
        subgroup_natural = new_subgroup(taxid, score)
        subgroup_artificial = new_subgroup(None, None)
    # Return Group object
    return Group(
        group_id=group_id,
        all=subgroup_all,
        natural=subgroup_natural,
        artificial=subgroup_artificial,
    )

def update_group(
        group: Group,
        taxid: int,
        score: float,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        artificial_taxids: set[int],
        ) -> tuple[Group, dict[int, list[int]]]:
    """
    Update a Group object with new taxid and score information.
    Args:
        group (Group): Group object to update.
        taxid (int): Taxid of new entry.
        score (float): Score of new entry.
        child_to_parent (dict[int, int]): Dictionary mapping each taxid to its parent.
        path_cache (dict[int, list[int]]): Cache of precomputed paths to the root.
        artificial_taxids (set[int]): Set of taxids that represent artificial or engineered sequences.
    Returns:
        tuple[Group, dict[int, list[int]]]: Tuple containing the updated
            Group object and the updated path cache.
    """
    # Determine if new entry is artificial
    is_artificial = taxid in artificial_taxids
    # Update all subgroup
    group.all, path_cache = update_subgroup(group.all, taxid, score, child_to_parent, path_cache)
    # Update natural and artificial subgroups as appropriate
    if is_artificial:
        group.artificial, path_cache = update_subgroup(group.artificial, taxid, score,
                                                       child_to_parent, path_cache)
    else:
        group.natural, path_cache = update_subgroup(group.natural, taxid, score,
                                        child_to_parent, path_cache)
    # Return updated Group object
    return group, path_cache

def output_group(
        group: Group,
        ) -> list[str]:
    """
    Process a Group object into a list of fields for output.
    """
    return [group.group_id] + output_subgroup(group.all) + \
        output_subgroup(group.natural) + output_subgroup(group.artificial)

#=======================================================================
# Functions for processing input and output
#=======================================================================

def get_output_header() -> list[str]:
    """
    Get the header for the output TSV.
    Returns:
        str: Header for the output TSV.
    """
    fields_base = ["lca", "n_entries", "top_taxid",
                   "min_score", "max_score", "mean_score"]
    fields_all = [f + "_all" for f in fields_base]
    fields_natural = [f + "_natural" for f in fields_base]
    fields_artificial = [f + "_artificial" for f in fields_base]
    return ["group_id"] + fields_all + fields_natural + fields_artificial

def write_output_line(
        fields: list[str],
        headers: list[str],
        output_file: TextIO,
        ) -> None:
    """
    Write a line to the output file.
    """
    assert len(fields) == len(headers), \
        f"Number of fields ({len(fields)}) does not match number of headers ({len(headers)})."
    line = "\t".join(fields) + "\n"
    output_file.write(line)

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
    assert group_field in header_fields, f"Group field not found in header: {group_field}"
    assert taxid_field in header_fields, f"Taxid field not found in header: {taxid_field}"
    assert score_field in header_fields, f"Score field not found in header: {score_field}"
    # Return indices
    group_idx = header_fields.index(group_field)
    taxid_idx = header_fields.index(taxid_field)
    score_idx = header_fields.index(score_field)
    return group_idx, taxid_idx, score_idx

def process_input_line(
        fields: list[str],
        group_idx: int,
        taxid_idx: int,
        score_idx: int,
        group_info: Group | None,
        child_to_parent: dict[int, int],
        path_cache: dict[int, list[int]],
        artificial_taxids: set[int],
        ) -> tuple[Group, dict[int, list[int]]]:
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
        artificial_taxids (set[int]): Set of taxids that represent artificial or
            engineered sequences.
    Returns:
        tuple[Group, dict[int, list[int]]]: Tuple containing the updated
            Group object containing information about the corresponding entry
            group and the updated path cache.
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
        # NB: new_group returns a Group object only
        return new_group(group_id, taxid, score, artificial_taxids), path_cache
    # Otherwise, update Group object
    # NB: update_group returns a tuple with the updated Group object and the
    # updated path cache
    return update_group(group_info, taxid, score, child_to_parent,
                        path_cache, artificial_taxids)

def parse_input_tsv(
        input_path: str,
        output_path: str,
        group_field: str,
        taxid_field: str,
        score_field: str,
        child_to_parent: dict[int, int],
        artificial_taxids: set[int],
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
        artificial_taxids (set[int]): Set of taxids that represent artificial or
            engineered sequences.
    """
    with open_by_suffix(input_path) as inf, open_by_suffix(output_path, "w") as outf:
        # Write header to output file
        output_header = get_output_header()
        write_output_line(output_header, output_header, outf)
        # Read header from input file
        header = inf.readline().strip().split("\t")
        # Get indices of fields
        group_idx, taxid_idx, score_idx = parse_input_header(
            header, group_field, taxid_field, score_field)
        path_cache = {}
        def parse_line(fields: list[str], group_info: Group | None) -> tuple[Group, dict[int, list[int]]]:
            return process_input_line(fields, group_idx, taxid_idx, score_idx, group_info,
                                      child_to_parent, path_cache, artificial_taxids)
        # Process first line
        fields = inf.readline().strip().split("\t")
        group_id = fields[group_idx]
        group_info, path_cache = parse_line(fields, None)
        # Iterate over input file
        n_entries = 1
        n_groups = 1
        for line in inf:
            # Parse fields
            fields = line.strip().split("\t")
            # Get group ID and check sorting
            group_id_new = fields[group_idx]
            assert group_id_new >= group_id, \
                f"Group ID out of order: {group_id_new} < {group_id}"
            # If group ID is unchanged, process line into existing Group object
            if group_id_new == group_id:
                group_info, path_cache = parse_line(fields, group_info)
            # Otherwise, write previous Group object to output and start new one
            else:
                out_line = output_group(group_info)
                write_output_line(out_line, output_header, outf)
                group_info, path_cache = parse_line(fields, None)
                group_id = group_id_new
                n_groups += 1
            n_entries += 1
        # Write last Group object to output
        out_line = output_group(group_info)
        write_output_line(out_line, output_header, outf)
        logger.info(f"Processed {n_entries} entries in {n_groups} groups.")

#=======================================================================
# Miscellaneous functions
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
    parser.add_argument("--artificial", "-a",
                        help="Parent taxid for artificial sequences (to be handled separately).")
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

def parse_taxonomy_db(path: str, artificial_taxid: int
                      ) -> tuple[dict[int, int], dict[int, set[int]]]:
    """
    Parse taxonomy DB into two dictionaries: one mapping each taxid
    to its parent taxid, and one mapping each taxid to its children taxids.
    Args:
        path (str): Path to taxonomy DB.
        artificial_taxid (int): Taxid of artificial parent.
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
            taxid = int(fields[0])
            parent_taxid = int(fields[2])
            child_to_parent[taxid] = parent_taxid
            parent_to_children[parent_taxid].add(taxid)
    # Check that DB contains root as the topmost taxid
    assert TAXID_ROOT in child_to_parent and TAXID_ROOT in parent_to_children, \
        "Taxonomy DB does not contain root."
    assert child_to_parent[TAXID_ROOT] == TAXID_ROOT, "Root taxid has a parent."
    assert TAXID_ROOT in parent_to_children[TAXID_ROOT], "Root taxid must be its own child."
    # Check that DB contains artificial parent taxid
    assert artificial_taxid in child_to_parent, \
        f"Artificial parent taxid not found in child-to-parent DB: {artificial_taxid}, {child_to_parent}"
    assert artificial_taxid in parent_to_children, \
        f"Artificial parent taxid not found in parent-to-children DB: {artificial_taxid}, {parent_to_children}"
    # Return dictionaries
    return child_to_parent, parent_to_children

def get_descendants(taxid: int, parent_to_children: dict[int, set[int]]) -> set[int]:
    """
    Get a set of all descendants of a taxid.
    Args:
        taxid (int): Taxid to get descendants of.
        parent_to_children (dict[int, set[int]]): Dictionary mapping each taxid
            to its children taxids.
    Returns:
        set[int]: Set of all descendants of the taxid, including the taxid itself.
    """
    descendants = set([taxid])
    descendants_new = parent_to_children[taxid]
    while descendants_new:
        descendants.update(descendants_new)
        descendants_newer = set()
        for taxid in descendants_new:
            descendants_newer.update(parent_to_children[taxid])
        # Remove any children already in set
        descendants_newer = descendants_newer - descendants
        # Update set of children
        descendants_new = descendants_newer
    # Return set of descendants
    return descendants

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
    child_to_parent, parent_to_children = parse_taxonomy_db(args.db, int(args.artificial))
    logger.info(f"Parsed taxonomy information for {len(child_to_parent)} taxids.")
    # Get set of artificial taxids
    logger.info("Getting set of artificial taxids.")
    artificial_taxids = get_descendants(int(args.artificial), parent_to_children)
    logger.info(f"Found {len(artificial_taxids)} artificial taxids.")
    # Parse input TSV and write LCA information to output TSV
    logger.info("Parsing input TSV.")
    parse_input_tsv(args.input, args.output, args.group, args.taxid, args.score,
                    child_to_parent, artificial_taxids)
    # Log completion
    logger.info("Script complete.")

if __name__ == "__main__":
    main()