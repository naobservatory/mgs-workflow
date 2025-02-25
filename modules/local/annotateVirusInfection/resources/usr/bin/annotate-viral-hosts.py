#!/usr/bin/env python

#=======================================================================
# Preamble
#=======================================================================

# Import libraries
import pandas as pd
import subprocess
from typing import Dict, Set, List
import argparse
from collections import defaultdict
import logging
from datetime import datetime, timezone

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
# Auxiliary functions
#=======================================================================

def get_virus_host_mapping(db_path: str) -> Dict[str, Set[str]]:
    """
    Import a TSV from Virus-Host-DB and extract virus/host info into
    a dictionary.

    Args:
        db_path (str): Path to the Virus-Host-DB TSV File.

    Returns:
        Dict[str, Set[str]]: A dictionary mapping virus taxids to
            sets of host taxids.
    """
    logger.info("Importing Virus-Host-DB.")
    df = pd.read_csv(db_path, sep="\t", dtype=str)
    df_cleaned = df.loc[df["host tax id"].notnull()]
    logger.info("Generating mapping from Virus-Host-DB.")
    mapping = df.groupby("virus tax id")["host tax id"].apply(set).to_dict()
    return mapping

def expand_taxid(taxid: str, nodes: pd.DataFrame) -> Set[str]:
    """
    Given a starting taxid and an NCBI nodes DataFrame, return a set of
    all taxids descended from that starting taxid (including itself).

    Args:
        taxid (str): Ancestral taxid as a numeric string.
        nodes (pd.DataFrame): A DataFrame generated from an NCBI nodes DMP
            file giving the parent taxid of each taxid within the NCBI

    Returns:
        Set[str]: Set of descendant taxids, including the original
            ancestor.
    """
    taxids_old = {taxid}
    taxids_desc = nodes.loc[nodes["parent_taxid"].isin(taxids_old)]["taxid"]
    taxids_new = taxids_old.union(taxids_desc)
    while taxids_new > taxids_old:
        taxids_old = taxids_new
        taxids_desc = nodes.loc[nodes["parent_taxid"].isin(taxids_old)]["taxid"]
        taxids_new = taxids_old.union(taxids_desc)
    return taxids_new

def get_host_taxids(hosts: Dict[str, str], nodes: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Given a dictionary mapping host group names to host ancestor taxids,
    return a dictionary instead mapping those host names to sets of
    descendant taxids.

    Args:
        hosts (Dict[str, str]): A dictionary mapping string names for
            host taxa (e.g. "human", "vertebrate") to ancestral taxids
            as numeric strings (e.g. "9606", "7742").
        nodes (pd.DataFrame): A DataFrame generated from an NCBI nodes DMP
            file giving the parent taxid of each taxid within the NCBI
            taxonomy structure.

    Returns:
        Dict[str, Set[str]]: A dictionary mapping those same string
            names to sets of descendant taxids returned by
            expand_taxid().
    """
    logger.info("Generating host taxid dictionary from input TSV.")
    host_dict_out = dict()
    for k in hosts.keys():
        host_dict_out[k] = expand_taxid(hosts[k], nodes)
        n_desc = len(host_dict_out[k])
        logger.info(f"\tFound {n_desc} descendant taxids for host group \"{k}\" (taxid {hosts[k]}).")
    logger.info("Host taxid dictionary construction complete.")
    return host_dict_out

def build_virus_tree(viral_taxa_df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Build a dictionary representing the viral taxonomic tree from a
    viral taxa DataFrame.

    Args:
        viral_taxa_df (pd.DataFrame): DataFrame containing viral
            taxa information. Must have taxid and parent_taxid
            columns.

    Returns:
        Dict[str, Set[str]]: A dictionary mapping virus taxids 
            to sets of their child taxids.
    """
    logger.info("Generating viral parent-child tree from viral taxa DB.")
    tree = defaultdict(set)
    for _, row in viral_taxa_df.iterrows():
        tree[row['parent_taxid']].add(row['taxid'])
    logger.info("Viral parent-child tree generation complete.")
    return tree

def mark_direct_infections(virus_taxids: pd.Series,
                           host_taxids: Set[str],
                           virus_host_mapping: Dict[str, Set[str]]) -> pd.Series:
    """For a given set of host taxids, check whether each virus taxid in a
    series is directly marked as infecting any of those host taxids,
    according to Virus-Host-DB. Doesn't consider the infection status of
    ancestral or descendant taxids.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        host_taxids (Set[str]): Series of host taxids to check against.
            Should be expanded to include all descendant taxids using
            netch_descendants().
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().

    Returns:
        pd.Series: Integer series giving direct infection status of each
            viral taxid in virus_taxids.
    """
    def check_direct_infection(virus_taxid):
        if virus_taxid in virus_host_mapping:
            return int(bool(virus_host_mapping[virus_taxid] & host_taxids))
        return -1 # Indicates not in mapping; needs special handling downstream
    statuses = virus_taxids.apply(check_direct_infection)
    return pd.Series(statuses.values, index=virus_taxids)

def add_descendants(virus_tree: Dict[str, Set[str]],
                    taxids_start: Set[str]) -> Set[str]:
    """
    Given a viral taxonomy tree and a set of starting taxids, return a set
    containing those taxids and all their descendants.

    Args:
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        taxids_start (Set[str]): Initial set of taxa for which to find
            descendants.

    Returns:
        Set[str]: Expanded set of taxids including taxids_start and all
            their descendants.
    """
    def expand_taxids(taxids):
        new_taxids = set(taxids)
        for taxid in taxids:
            new_taxids.update(virus_tree.get(taxid, set()))
        return new_taxids
    taxids = taxids_start
    taxids_new = expand_taxids(taxids)
    while taxids_new > taxids:
        taxids = taxids_new
        taxids_new = expand_taxids(taxids)
    return taxids_new

def mark_descendant_infections(virus_tree: Dict[str, Set[str]],
                               statuses: pd.Series) -> pd.Series:
    """
    Given a series of infection statuses, propagate statuses
    from each taxid to its descendants and return a new status series.

    Args:
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        statuses (pd.Series): Integer series of infection statuses
            for virus_taxids.

    Returns:
        pd.Series: Integer series of updated infection statuses.
    """
    # All descendents of a taxid marked 1 should be marked 1
    taxids_1 = set(statuses.index[statuses == 1])
    taxids_1_expanded = add_descendants(virus_tree, taxids_1)
    statuses[statuses.index.isin(taxids_1_expanded)] = 1
    # All descendents of a taxid marked 0 should be marked 0
    taxids_0 = set(statuses.index[statuses == 0])
    taxids_0_expanded = add_descendants(virus_tree, taxids_0)
    statuses[statuses.index.isin(taxids_0_expanded)] = 0
    # Descendants of a taxid marked 3 should be marked 3 iff they are not
    # already marked 1
    taxids_3 = set(statuses.index[statuses == 3])
    taxids_3_expanded = add_descendants(virus_tree, taxids_3)
    assert statuses[statuses.index.isin(taxids_3_expanded)].isin([0, 2]).sum() == 0, \
            "Some taxids marked 3 have descendants marked 0 or 2."
    statuses[statuses.index.isin(taxids_3_expanded) & ~statuses.isin([1])] = 3
    # Descendants of a taxid marked 2 should be marked 2 iff they are marked -1
    taxids_2 = set(statuses.index[statuses == 2])
    taxids_2_expanded = add_descendants(virus_tree, taxids_2)
    statuses[statuses.index.isin(taxids_2_expanded) & statuses.isin([-1])] = 2
    # Finally, mark any remaining -1 taxids as 0
    statuses[statuses == -1] = 0
    # After this, there should be no taxids marked -1
    assert statuses.isin([-1]).sum() == 0, "Some taxids are still unresolved (marked -1)."
    return statuses

def exclude_infections(virus_tree: Dict[str, Set[str]],
                       statuses: pd.Series,
                       exclude_taxids: List[str]) -> pd.Series:
    """
    Given a taxonomy tree, a series of infection status generated by
    mark_descendant_infections(), and a set of taxids to exclude,
    assign all taxids descended from those taxids a status of 0.

    Args:
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        statuses (pd.Series): Integer series of infection statuses
            for virus_taxids.
        exclude_taxids (List[str]): List of virus taxids to exclude from
            host annotation (i.e. to force to mark as non-infecting).

    Returns:
        pd.Series: Integer series of infection statuses, updated to
            exclude the specified taxa by marking them 0.
    """
    exclude_taxids_expanded = add_descendants(virus_tree, set(exclude_taxids))
    statuses[statuses.index.isin(exclude_taxids_expanded)] = 0
    return statuses

def mark_ancestor_infections_single(virus_taxid: str,
                                    virus_df: pd.DataFrame,
                                    virus_tree: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Given a single viral taxid of interest, a DataFrame giving current
    infection statuses for each of a set of viral taxids, and
    a mapping of each taxid to its children, conduct a recrusive, 
    memoized tree search to determine the infection of the focal taxid
    based on the infection statuses of its descendants.

    Args:
        virus_taxid (str): Specific viral taxid whose status to check.
        virus_df (pd.DataFrame): DataFrame with a column "taxid" giving
            a list of all viral taxids, a column "status" giving
            their currently-assigned infection status (0, 1 or 2),
            and a boolean column "checked" indicating whether its status
            has been checked and assigned by this function.
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree
            generated from a viral DB using build_virus_tree(),
            mapping each taxid to a set of its child taxids.

    Returns:
        pd.DataFrame: An updated DataFrame in the same format as 
            virus_df.
    """
    # Define auxiliary functions
    def set_checked(taxid, df):
        df.loc[taxid,"checked"] = True
        return df
    def set_status(taxid, df, value):
        df.loc[taxid,"status"] = value
        return set_checked(taxid, df)
    # Before passing to this function, already confirmed that taxid
    # has children and is unresolved, so start by getting children
    children = virus_df.loc[list(virus_tree[virus_taxid])]
    # If all children have been checked already, evaluate on the basis
    # of child statuses
    if children["checked"].all():
        child_statuses = children["status"]
        # If all children are marked -1, keep current status
        if (child_statuses == -1).all():
            return set_checked(virus_taxid, virus_df)
        # Otherwise, if all children have the same status, inherit that status
        elif (child_statuses == child_statuses.iloc[0]).all():
            return set_status(virus_taxid, virus_df, child_statuses.iloc[0])
        # If any child is marked 2, mark this taxid as 2
        elif child_statuses.isin([2]).any():
            return set_status(virus_taxid, virus_df, 2)
        # If there are both children marked 0 and marked [1 or 3], mark this taxid as 2
        elif (child_statuses.isin([0]).any() and
              child_statuses.isin([1,3]).any()):
            return set_status(virus_taxid, virus_df, 2)
        # If any child is marked 0, keep current status
        elif (child_statuses == 0).any():
            return set_checked(virus_taxid, virus_df)
        # Otherwise, mark as 3
        # This covers all mixtures of 1, -1 and 3
        else:
            return set_status(virus_taxid, virus_df, 3)
    # Otherwise, run this function for each unchecked child, then repeat
    children_unchecked = children.loc[~children["checked"]]
    for child in children_unchecked.index:
        virus_df = mark_ancestor_infections_single(child, virus_df, virus_tree)
    return mark_ancestor_infections_single(virus_taxid, virus_df, virus_tree)

def mark_ancestor_infections(virus_taxids: pd.Series,
                             virus_tree: Dict[str, Set[str]],
                             statuses: pd.Series) -> pd.Series:
    """
    Given a series of infection statuses, propagate statuses
    to ancestor taxa and return a new status series.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree
            generated from a viral DB using build_virus_tree(),
            mapping each taxid to a set of its child taxids.
        statuses (pd.Series): Integer series of direct infection statuses
            for virus_taxids.

    Returns:
        pd.Series: Integer series of updated infection statuses.
    """
    # Create a dataframe of taxids and current statuses
    df = pd.DataFrame({"status": statuses})
    # Mark taxids maked 1 as checked
    df["checked"] = df["status"] == 1
    # Get number of children and mark childless nodes as checked
    df["n_children"] = df.index.map(lambda x: len(virus_tree[x]))
    df["checked"] = df["checked"] | (df["n_children"] == 0)
    # For nodes with children, count how many have been checked
    # Iterate over unchecked rows until all rows have been checked
    while not df["checked"].all():
        for taxid in df.index:
            if df.loc[taxid]["checked"]:
                continue
            df = mark_ancestor_infections_single(taxid, df, virus_tree)
    return df["status"]

def check_infection(virus_taxids: pd.Series,
                    host_taxids: Set[str],
                    virus_tree: Dict[str, Set[str]],
                    virus_host_mapping: Dict[str, Set[str]],
                    hard_exclude_taxids: List[str]) -> pd.Series:
    """
    For a single set of host taxids, checks whether each virus taxid in a
    series infects that group of hosts and classifies that taxid according
    to its infection status:
        - 1 designates taxa that are affirmatively marked as infecting that
          host group. It includes taxa that are directly marked in Virus-Host
          DB as infecting that host group; taxa that are descended from such
          taxa; and taxa whose descendants are all marked 1.
        - 0 designates taxa that are marked as not infecting that host group.
          It includes taxa that: are included in Virus-Host DB; are not marked
          in Virus-Host DB as infecting that host group; and have no 1-marked
          descendants. It also includes taxa that are descended from such taxa,
          and taxa whose descendants are all marked 0.
        - 3 designates taxa that are likely to infect that host group, but cannot
          be affirmatively marked as such. It includes taxa that are not included
          in Virus-Host DB, and have no 0-marked descendants, but cannot be marked
          1 according to the rules above.
        - 2 designates taxa of uncertain infection status. It includes all other taxa:
          those that are not included in Virus-Host DB and have both 0- and 1-marked
          descendants, and those descended from a 2-marked taxon that are not given
          another status by the rules above.
    If a taxon is included in hard_exclude_taxids, it and all descendants
    are marked 0 regardless of Virus-Host-DB; this is intended to account
    for misannotated taxa.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        host_taxids (Set[str]): Series of host taxids to check against.
            Should be expanded to include all descendant taxids using
            expand_taxid().
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().
        hard_exclude_taxids (List[str]): List of virus taxids to hard-exclude
            from host annotation (i.e. to force to mark as non-infecting).

    Returns:
        pd.Series: Four-state integer series of infection statuses for each viral 
            taxid in virus_taxids.
    """
    # Start by marking direct infections
    logger.info("\tMarking direct infection status from Virus-Host-DB.")
    statuses = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
    # Exclude hard-excluded taxa
    logger.info("\tExclude hardcoded non-infecting taxids.")
    statuses = exclude_infections(virus_tree, statuses, hard_exclude_taxids)
    # Expand to ancestors
    logger.info("\tPropagating infection status to ancestor taxids.")
    statuses = mark_ancestor_infections(virus_taxids, virus_tree, statuses)
    # Expand to descendants
    logger.info("\tPropagating infection status to descendant taxids.")
    statuses = mark_descendant_infections(virus_tree, statuses)
    logger.info("\tInfection status inference complete.")
    return statuses

def annotate_virus_db_single(virus_db: pd.DataFrame,
                             host_name: str,
                             host_taxids: Set[str],
                             virus_tree: Dict[str, Set[str]],
                             virus_host_mapping: Dict[str, Set[str]],
                             hard_exclude_taxids: List[str]) -> pd.DataFrame:
    """
    Annotate a DataFrame of virus taxa with infection status information
    for a specified host taxon.

    Args:
        virus_db (pd.DataFrame): DataFrame with columns indicating taxid
            and parent taxid for each virus.
        host_name (str): String name for host taxon of interest, used to
            label infection status column.
        host_taxids (Set[str]): Expanded set of taxids within host taxon,
            generated with get_host_taxids().
        virus_tree (Dict[str, Set[str]]): Dictionary giving viral
            taxonomic structure, generated with build_virus_tree().
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().
        hard_exclude_taxids (List[str]): List of virus taxids to hard-exclude
            from host annotation (i.e. to force to mark as non-infecting).

    Returns:
        pd.DataFrame: Copy of virus_db with an additional column giving
            infection status for the specified host taxon.
    """

    # Get infection statuses
    logger.info(f"Inferring infection statuses for {host_name}-infecting viruses.")
    statuses = check_infection(virus_db["taxid"], host_taxids, virus_tree,
                               virus_host_mapping, hard_exclude_taxids)
    logger.info(f"Infection statuses for {host_name}-infecting viruses inferred.")
    # Annotate DB with infection statuses and return
    virus_db.set_index("taxid", inplace=True, drop=False)
    virus_db["infection_status_" + host_name] = statuses
    #virus_db.reset_index(inplace=True, drop=True)
    logger.info(f"Virus DB annotated with {host_name}-infection status.")
    return virus_db

def annotate_virus_db(virus_db: pd.DataFrame,
                      host_mapping: Dict[str, Set[str]],
                      virus_host_mapping: Dict[str, Set[str]],
                      hard_exclude_taxids: List[str]) -> pd.DataFrame:
    """
    Given a DataFrame of virus taxa (including taxids and parent taxids)
    and another of host taxa (including names and taxids) add a column
    to the former for each taxid in the latter indicating infection
    status for each virus.

    Args:
        virus_db (pd.DataFrame): DataFrame with columns indicating taxid
            and parent taxid for each virus.
        host_mapping (Dict[str, Set[str]]): Mapping of host names (e.g.
            "human", "bird") to sets of member taxids.
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().
        hard_exclude_taxids (List[str]): List of virus taxids to hard-exclude
            from host annotation (i.e. to force to mark as non-infecting).

    Returns:
        pd.DataFrame: Copy of virus_db with additional columns giving
            infection status for each taxon in host_db.
    """
    # Get viral taxonomic tree
    virus_tree = build_virus_tree(virus_db)
    # Add annotations for each host group
    for k in host_mapping.keys():
        virus_db = annotate_virus_db_single(virus_db, k, host_mapping[k],
                                            virus_tree, virus_host_mapping,
                                            hard_exclude_taxids)
    return virus_db

#=======================================================================
# Main function
#=======================================================================

def main():
    logger.info("Initializing script.")
    # Define argument parsing
    parser = argparse.ArgumentParser(description="Annotate a viral TSV with host infection status.")
    parser.add_argument("virus_db", help="Path to TSV of virus taxids and tree structure.")
    parser.add_argument("host_db", help="Path to TSV of target host taxids and names.")
    parser.add_argument("infection_db", help="Path to TSV from Virus-Host DB giving host information.")
    parser.add_argument("nodes_db", help="Path to NCBI taxonomy nodes file.")
    parser.add_argument("hard_exclude_taxids", help="Space-delimited list of viral taxids to hard-exclude from host annotation.")
    parser.add_argument("output_db", help="Output path for host-annotated TSV.")
    args = parser.parse_args()
    # Import inputs
    logger.info("Importing input TSVs.")
    virus_db = pd.read_csv(args.virus_db, sep="\t", dtype=str)
    host_db = pd.read_csv(args.host_db, sep="\t", dtype=str)
    nodes_db = pd.read_csv(args.nodes_db, sep="\t", dtype=str, header=None
                           ).iloc[:,[0,2]].rename(columns={0:"taxid",2:"parent_taxid"})
    virus_host_mapping = get_virus_host_mapping(args.infection_db)
    hard_exclude_taxids = args.hard_exclude_taxids.split(" ")
    # Prepare dictionary of host taxids
    host_dict_single = host_db.set_index("name")["taxid"].to_dict()
    host_dict_full = get_host_taxids(host_dict_single, nodes_db)
    # Add annotations
    logger.info("Initializing infection-state annotation.")
    output_db = annotate_virus_db(virus_db, host_dict_full, virus_host_mapping, 
                                  hard_exclude_taxids)
    # Write output
    logger.info("Writing output.")
    output_db.to_csv(args.output_db, sep="\t", index=False)
    logger.info("Script complete.")

if __name__ == "__main__":
    main()
