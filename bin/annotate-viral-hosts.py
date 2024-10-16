#!/usr/bin/env python

#=======================================================================
# Preamble
#=======================================================================

# Import libraries
import pandas as pd
import subprocess
from typing import Dict, Set
import argparse
from collections import defaultdict

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
    df = pd.read_csv(db_path, sep="\t")
    mapping = df.groupby("virus tax id")["host tax id"].apply(set).to_dict()
    return mapping

def run_gimme_taxa(taxid: str) -> str:
    """
    Run gimme_taxa.py for a given taxid and capture stderr output.

    Args:
        taxid (str): Target taxid as a numeric string.

    Returns:
        str: TSV output as a tab- and newline-delimited string.
    """
    cmd = ["gimme_taxa.py", str(taxid)]
    try:
        p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return p.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching descendants for taxid {taxid}: {e}")
        print(f"Error output: {e.stderr}")
        raise  # Re-raise the exception if it's not about missing database

def fetch_descendents(taxid: str) -> Set[str]:
    """
    For a given taxid, run gimme_taxa.py and process the output into
    list format.

    Args:
        taxid (str): Target taxid as a numeric string.

    Returns:
        Set[str]: Set of descendent taxids (including the original
            taxid as the first element).
    """
    # Run gimme_taxa.py
    try:
        output = run_gimme_taxa(taxid)
    except Exception as e:
        print(f"Error: {e}")
        raise  # Re-raise the exception to be handled by the calling function
    desc_split = output.split("\n")
    # Get list of descendents
    descendents = [taxid]
    tab = False # Skip frontmatter from stdout
    for line in desc_split:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("parent_taxid"):
            tab = True # Stop skipping after header line
            continue
        elif not tab:
            continue
        try:
            parent_taxid, descendent_taxid, descendent_name = line.split("\t")
            descendents.append(int(descendent_taxid))
        except ValueError:
            print(line)
            raise
    return set(descendents)

def get_host_taxids(hosts: Dict[str, str]) -> Dict[str, Set[str]]:
    """
    Given a dictionary mapping host group names to host ancestor taxids,
    return a dictionary instead mapping those host names to sets of
    descendent taxids.

    Args:
        hosts (Dict[str, str]): A dictionary mapping string names for
            host taxa (e.g. "human", "vertebrate") to ancestral taxids
            as numeric strings (e.g. "9606", "7742").

    Returns:
        Dict[str, Set[str]]: A dictionary mapping those same string
            names to sets of descendent taxids returned by
            fetch_descendents().
    """
    host_dict_out = dict()
    for k in hosts.keys():
        host_dict_out[k] = fetch_descendents(hosts[k])
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
    tree = defaultdict(set)
    for _, row in viral_taxa_df.iterrows():
        tree[row['parent_taxid']].add(row['taxid'])
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
            Should be expanded to include all descendent taxids using
            fetch_descendents().
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().

    Returns:
        pd.Series: Binary series giving direct infection status of each
            viral taxid in virus_taxids.
    """
    def check_direct_infection(virus_taxid):
        if virus_taxid in virus_host_mapping:
            return int(bool(virus_host_mapping[virus_taxid] & host_taxids))
        return 0
    return virus_taxids.apply(check_direct_infection)

def mark_descendant_infections(virus_taxids: pd.Series,
                               virus_tree: Dict[str, Set[str]],
                               statuses: pd.Series) -> pd.Series:
    """
    Given a binary series of direct infection status generated by
    mark_direct_infections(), mark any virus taxids descended from marked
    taxids as also infecting the same host group.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        statuses (pd.Series): Binary series of direct infection statuses
            for virus_taxids, generated by mark_direct_infections().

    Returns:
        pd.Series: Binary series corresponding to each taxid in
            virus_taxids, in which each taxid is marked 1 if (i) it was
            marked 1 in the input statuses, or (ii) it is descended
            from a taxid that was so marked.
    """
    # Get the set of initially marked taxids
    marked_taxids = set(virus_taxids[statuses == 1])
    # Expand to include all their descendants
    def expand_taxids(taxids):
        new_taxids = set(taxids)
        for taxid in taxids:
            new_taxids.update(virus_tree.get(taxid, set()))
        return new_taxids
    marked_taxids_new = expand_taxids(marked_taxids)
    while marked_taxids_new > marked_taxids:
        marked_taxids = marked_taxids_new
        marked_taxids_new = expand_taxids(marked_taxids)
    # Mark those as 1
    final_status = pd.Series(index=virus_taxids, data=0)
    final_status[final_status.index.isin(marked_taxids)] = 1
    return final_status

def mark_ancestor_infections_single(virus_taxid: str,
                                    virus_df: pd.DataFrame,
                                    virus_tree: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Given a single viral taxid of interest, a DataFrame giving current
    infection statuses for each of a set of viral taxids, and
    a mapping of each taxid to its children, conduct a memoized tree
    search to determine the infection status of the focal taxid as
    follows:
        - If the taxid is already marked with a nonzero value, it
          remains so marked.
        - If all of the taxid's descendants are marked 1, it is also
          marked 1.
        - Otherwise, if any of its descendants are marked 1 or 2, it is
          marked 2.
        - Otherwise, it is marked 0.
    Determining the status of a given taxon thus requires determining
    the status of all its descendants. To make future iterations of
    this process more efficient, the function also updates the status 
    of those taxids in the DataFrame as it goes.

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
    # If virus_taxid has already been checked, no change
    def check_checked(taxid, df):
        return df.loc[df["taxid"] == taxid, "checked"].iloc[0]
    if check_checked(virus_taxid, virus_df):
        return virus_df
    # If virus_taxid already has a nonzero value, marked as checked with
    # no other change
    def get_status(taxid, df):
        return df.loc[df["taxid"] == taxid, "status"].iloc[0]
    def set_checked(taxid, df):
        df.loc[df["taxid"] == taxid, "checked"] = True
        return df
    if get_status(virus_taxid, virus_df) != 0:
        return set_checked(virus_taxid, virus_df)
    # Otherwise, get a list of all child taxids
    children = virus_tree.get(virus_taxid, set())
    # If there are no children, mark as 0 and return
    def set_status(taxid, df, value):
        df.loc[df["taxid"] == taxid, "status"] = value
        return set_checked(taxid, df)
    if not children:
        return set_status(virus_taxid, virus_df, 0)
    # Otherwise, check if all children have been checked
    child_checked = [check_checked(child, virus_df) for child in children]
    # If so, evaluate on the basis of child statuses
    if all (child_checked):
        child_statuses = [get_status(child, virus_df) for child in children]
        if all(status == 1 for status in child_statuses):
            return set_status(virus_taxid, virus_df, 1)
        elif any(status in [1,2] for status in child_statuses):
            return set_status(virus_taxid, virus_df, 2)
        else:
            return set_status(virus_taxid, virus_df, 0)
    # Otherwise, run this function for each child and repeat
    for child in children:
        virus_df = mark_ancestor_infections_single(child, virus_df, virus_tree)
    return mark_ancestor_infections_single(virus_taxid, virus_df, virus_tree)

def mark_ancestor_infections(virus_taxids: pd.Series,
                             virus_tree: Dict[str, Set[str]],
                             statuses: pd.Series) -> pd.Series:
    """
    Given a binary series of direct & descendant infection status
    generated by mark_descendant_infections(), return a three-state
    series marking each virus taxid as follows:
        - If a taxid is already marked 1, it remains so marked.
        - If all of the taxon's descendent taxa are marked 1, it is also
          marked 1.
        - If some but not all of the taxon's descendents are marked 1,
          or if any is marked 2, it is marked 2.
        - Otherwise, the taxon is marked 0.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree
            generated from a viral DB using build_virus_tree(),
            mapping each taxid to a set of its child taxids.
        statuses (pd.Series): Binary series of direct infection statuses
            for virus_taxids, generated by mark_descendant_infections().

    Returns:
        pd.Series: Three-state series corresponding to each taxid in
            virus_taxids, marked as described above.
    """
    # Create a dataframe of taxids and current statuses
    df = pd.DataFrame({"taxid": virus_taxids, "status": statuses,
                       "checked": False})
    df.set_index("taxid", inplace=True)
    # Call mark_ancestor_infections_single() on each row
    for taxid in virus_taxids:
        df = mark_ancestor_infections_single(taxid, df, virus_tree)
    return df["status"]

def check_infection(virus_taxids: pd.Series,
                    host_taxids: Set[str],
                    virus_tree: Dict[str, Set[str]],
                    virus_host_mapping: Dict[str, Set[str]]) -> pd.Series:
    """
    For a single set of host taxids, check whether each virus taxid in a
    series infects that group of hosts. Classifies each virus taxon as
    follows:
        - If the taxon is directly marked in Virus-Host-DB as infecting
          that host group, returns 1.
        - If the taxon is descended from a taxon marked 1, it is also
          marked 1.
        - If all of the taxon's descendent taxa are marked 1, it is also
          marked 1.
        - If some but not all of the taxon's descendents are marked 1,
          or if any are marked 2, it is marked 2.
        - Otherwise, the taxon is marked 0.

    Args:
        virus_taxids (pd.Series): Series of virus taxids to check.
        host_taxids (Set[str]): Series of host taxids to check against.
            Should be expanded to include all descendent taxids using
            fetch_descendents().
        virus_tree (Dict[str, Set[str]]): Viral taxonomic tree,
            generated from a viral DB using build_virus_tree().
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().

    Returns:
        pd.Series: Series of infection statuses for each viral taxid in
            virus_taxids.
    """
    # Start by marking direct infections
    statuses = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
    # Expand to descendants
    statuses = mark_descendant_infections(virus_taxids, virus_tree, statuses)
    # Expand to ancestors and return
    statuses = mark_ancestor_infections(virus_taxids, virus_tree, statuses)
    return statuses

def annotate_virus_db_single(virus_db: pd.DataFrame,
                             host_name: str,
                             host_taxids: Set[str],
                             virus_tree: Dict[str, Set[str]],
                             virus_host_mapping: Dict[str, Set[str]]) -> pd.DataFrame:
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

    Returns:
        pd.DataFrame: Copy of virus_db with an additional column giving
            infection status for the specified host taxon.
    """

    # Get infection statuses
    statuses = check_infection(virus_db["taxid"], host_taxids, virus_tree, 
                               virus_host_mapping)
    # Annotate DB with infection statuses and return
    virus_db["infection_status_" + host_name] = statuses
    return virus_db

def annotate_virus_db(virus_db: pd.DataFrame,
                      host_db: pd.DataFrame,
                      virus_host_mapping: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Given a DataFrame of virus taxa (including taxids and parent taxids)
    and another of host taxa (including names and taxids) add a column
    to the former for each taxid in the latter indicating infection
    status for each virus.

    Args:
        virus_db (pd.DataFrame): DataFrame with columns indicating taxid
            and parent taxid for each virus.
        host_db (pd.DataFrame): DataFrame with columns indicating taxid
            and name for each host taxon of interest.
        virus_host_mapping (Dict[str, Set[str]]): Mapping of virus
            taxids to host taxids, generated using
            get_virus_host_mapping().

    Returns:
        pd.DataFrame: Copy of virus_db with additional columns giving
            infection status for each taxon in host_db.
    """
    # Get viral taxonomic tree
    virus_tree = build_virus_tree(virus_db)
    # Get dictionary of host taxid sets
    host_dict_single = host_db.set_index("name")["taxid"].to_dict()
    host_dict_full = get_host_taxids(host_dict_single)
    # Add annotations for each host group
    for k in host_dict_full.keys():
        virus_db = annotate_virus_db_single(virus_db, k, host_dict_full[k],
                                            virus_tree, virus_host_mapping)
    return virus_db

#=======================================================================
# Main function
#=======================================================================

def main():
    # Define argument parsing
    parser = argparse.ArgumentParser(description="Annotate a viral TSV with host infection status.")
    parser.add_argument("virus_db", help="Path to TSV of virus taxids and tree structure.")
    parser.add_argument("host_db", help="Path to TSV of target host taxids and names.")
    parser.add_argument("infection_db", help="Path to TSV from Virus-Host DB giving host information.")
    parser.add_argument("output_db", help="Output path for host-annotated TSV.")
    args = parser.parse_args()
    # Import inputs
    virus_db = pd.read_csv(args.virus_db, sep="\t")
    host_db = pd.read_csv(args.host_db, sep="\t")
    virus_host_mapping = get_virus_host_mapping(args.infection_db)
    # Add annotations
    output_db = annotate_virus_db(virus_db, host_db, virus_host_mapping)
    # Write output
    output_db.to_csv(args.output_db, sep="\t")

if __name__ == "__main__":
    main()
