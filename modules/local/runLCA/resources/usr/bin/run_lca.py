#!/usr/bin/env python

import sys
from typing import List, Dict, Optional, Any

# Define the root taxid for the NCBI taxonomy
NCBI_ROOT_TAXID = 1


def load_parent_map(tree_filepath: str) -> Dict[int, int]:
    """
    Loads the taxonomy parent map from a file.

    Assumes the file is tab-separated with taxid and parent_taxid.
    Handles the root node where parent_taxid is often the same as taxid.

    Args:
        tree_filepath: Path to the taxonomy tree file.

    Returns:
        A dictionary mapping taxid (int) to parent_taxid (int).
        Returns an empty dictionary if the file cannot be read.
    """
    parent_map: Dict[int, int] = {}
    try:
        with open(tree_filepath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        taxid = int(parts[0])
                        parent_taxid = int(parts[2])
                        parent_map[taxid] = parent_taxid
                    except ValueError:
                        sys.stderr.write(
                            f"Warning: Skipping malformed line in tree file: {line.strip()}\n"
                        )
                        pass
    except FileNotFoundError:
        sys.stderr.write(f"Error: Tree file not found at {tree_filepath}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"Error loading tree file: {e}\n")
        return {}

    if NCBI_ROOT_TAXID not in parent_map:
        parent_map[NCBI_ROOT_TAXID] = NCBI_ROOT_TAXID  # Parent of root is root

    return parent_map


def load_seq_taxids(
    seq_taxid_filepath: str, group_col_name: str, taxid_col_name: str
) -> Dict[Any, List[int]]:
    """
    Loads the mapping from seq_id to a list of taxids.

    Assumes the file is tab-separated with seq_id and taxid.

    Args:
        seq_taxid_filepath: Path to the seq_id/taxid mapping file.

    Returns:
        A dictionary mapping seq_id (can be string or int) to a list of taxids (int).
        Returns an empty dictionary if the file cannot be read.
    """
    seq_taxids: Dict[Any, List[int]] = {}
    try:
        with open(seq_taxid_filepath, "r") as f:
            first_line = f.readline()
            if not first_line.strip():
                pass  # Empty file case
            else:
                header_cols = first_line.strip().lower().split("\t")
                if taxid_col_name.lower() not in header_cols:
                    raise ValueError(f"Column does not exist: {taxid_col_name}")
                elif group_col_name.lower() not in header_cols:
                    raise ValueError(f"Column does not exist: {group_col_name}")
                else:
                    taxid_index = header_cols.index(taxid_col_name.lower())
                    group_index = header_cols.index(group_col_name.lower())
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    seq_id_str = parts[group_index]
                    try:
                        seq_id: Any = int(seq_id_str)
                    except ValueError:
                        seq_id = seq_id_str

                    try:
                        taxid = int(parts[taxid_index])
                        seq_taxids.setdefault(seq_id, []).append(taxid)
                    except ValueError:
                        sys.stderr.write(
                            f"Warning: Skipping line with non-integer taxid in seq file: {line.strip()}\n"
                        )
                        pass
    except FileNotFoundError:
        sys.stderr.write(f"Error: Seq/Taxid file not found at {seq_taxid_filepath}\n")
        return {}
    except Exception as e:
        sys.stderr.write(f"Error loading seq/taxid file: {e}\n")
        return {}

    return seq_taxids


def get_path_to_root(taxid: int, parent_map: Dict[int, int]) -> List[int]:
    """
    Gets the path from a given taxid up to the root (taxid 1).

    Args:
        taxid: The starting taxid (int).
        parent_map: The dictionary mapping taxid to parent_taxid.

    Returns:
        A list of taxids starting from the input taxid up to the root,
        inclusive. Returns an empty list if the taxid is not in the parent map.
        Example: [12345, 1234, 123, 1]
    """
    path: List[int] = []
    current_taxid = taxid

    # Traverse up until we reach the root or an unknown taxid
    while current_taxid in parent_map:
        path.append(current_taxid)
        if (
            current_taxid == parent_map[current_taxid]
        ):  # Check for root (parent points to self)
            break
        current_taxid = parent_map[current_taxid]

    return path


def find_lca_two(taxid1: int, taxid2: int, parent_map: Dict[int, int]) -> Optional[int]:
    """
    Finds the Lowest Common Ancestor (LCA) of two taxids.

    Args:
        taxid1: The first taxid (int).
        taxid2: The second taxid (int).
        parent_map: The dictionary mapping taxid to parent_taxid.

    Returns:
        The LCA taxid (int), or None if inputs are invalid
        (e.g., not found in parent_map or paths don't connect to root,
        although the latter is unlikely with a valid tree fragment).
    """
    if taxid1 == taxid2:
        return taxid1

    path1 = get_path_to_root(taxid1, parent_map)
    path2 = get_path_to_root(taxid2, parent_map)

    # If either path is empty, one of the taxids was invalid/not in map
    if not path1 or not path2:
        sys.stderr.write(
            f"Warning: Could not get path for taxid {taxid1} or {taxid2}. Returning None for LCA.\n"
        )
        return None

    # Reverse paths to start from the root
    path1.reverse()
    path2.reverse()

    lca: Optional[int] = None

    # Walk down the paths until they diverge
    for n1, n2 in zip(path1, path2):
        if n1 == n2:
            lca = n1  # Keep track of the last common node
        else:
            break  # Paths diverge

    return lca


def find_lca_set(taxid_list: List[int], parent_map: Dict[int, int]) -> Optional[int]:
    """
    Finds the Lowest Common Ancestor (LCA) for a list of taxids.

    Uses iterative binary LCA calculation: LCA(t1, t2, t3) = LCA(LCA(t1, t2), t3).

    Args:
        taxid_list: A list of taxids (int).
        parent_map: The dictionary mapping taxid to parent_taxid.

    Returns:
        The LCA taxid (int) for the set, or None if the list is empty
        or if an LCA cannot be determined (e.g., due to invalid taxids).
    """
    if not taxid_list:
        return None
    if len(taxid_list) == 1:
        return taxid_list[0]

    # Start with the LCA of the first two taxids
    current_lca = find_lca_two(taxid_list[0], taxid_list[1], parent_map)

    if current_lca is None:
        sys.stderr.write(
            f"Warning: Initial LCA calculation failed for {taxid_list[:2]}. Returning None for the set.\n"
        )
        return None

    for next_taxid in taxid_list[2:]:
        current_lca = find_lca_two(current_lca, next_taxid, parent_map)
        if current_lca is None:
            sys.stderr.write(
                f"Warning: Subsequent LCA calculation failed with taxid {next_taxid}. Returning None for the set.\n"
            )
            return None

    return current_lca


def main(
    tree_filepath: str,
    seq_taxid_filepath: str,
    group_col_name: str,
    taxid_col_name: str,
    outfile_name: str,
):
    """
    Main function to load data, compute LCAs, and print results.

    Args:
        tree_filepath: Path to the taxonomy tree file.
        seq_taxid_filepath: Path to the seq_id/taxid mapping file.
    """
    print("Loading taxonomy tree...")
    parent_map = load_parent_map(tree_filepath)
    if not parent_map:
        sys.stderr.write("Failed to load taxonomy tree. Exiting.\n")
        sys.exit(1)
    print(f"Loaded parent map with {len(parent_map)} entries.")

    print("Loading seq_id to taxids mapping...")
    seq_taxids = load_seq_taxids(seq_taxid_filepath, group_col_name, taxid_col_name)
    if not seq_taxids:
        sys.stderr.write("Failed to load seq_id to taxids mapping. Exiting.\n")
        sys.exit(1)
    print(f"Loaded mappings for {len(seq_taxids)} seq_ids.")

    print("Computing LCAs...")

    with open(outfile_name, "w") as file:
        for seq_id, taxid_list in seq_taxids.items():
            if not taxid_list:
                file.write(f"{seq_id}\tNone\n")
                continue

            lca_taxid = find_lca_set(taxid_list, parent_map)

            file.write(f"{seq_id}\t{lca_taxid}\n")

    print("Processing complete. Results written to output.txt.")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(
            "Usage: python your_script_name.py <taxonomy_tree_file> <seq_taxid_file> <group_col_name> <taxid_col_name> <outfile_name>"
        )
        sys.exit(1)

    tree_file = sys.argv[1]
    seq_taxid_file = sys.argv[2]
    group_col_name = sys.argv[3]
    taxid_col_name = sys.argv[4]
    outfile_name = sys.argv[5]

    main(tree_file, seq_taxid_file, group_col_name, taxid_col_name, outfile_name)

