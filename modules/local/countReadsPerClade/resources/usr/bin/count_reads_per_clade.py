import argparse
import gzip
import bz2
import csv
from collections import Counter, defaultdict
from typing import Iterator, TextIO

TaxId = str
Tree = dict[TaxId, list[TaxId]]


def open_by_suffix(filename: str, mode: str = "r", debug: bool = False) -> TextIO:
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode + "t")
    elif filename.endswith(".bz2"):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


def read_tsv(file_path: str) -> Iterator[dict[str, str]]:
    """
    Read a TSV file and yield rows one at a time.

    Args:
        file_path (str): Path to the TSV file

    Yields:
        dict: Dictionary representing each row
    """
    with open_by_suffix(file_path, mode="rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            yield row


def is_duplicate(read: dict[str, str]) -> bool:
    """
    Check if a read is a duplicate based on sequence ID and primary alignment exemplar.

    Args:
        read: Dictionary representing a read record

    Returns:
        True if the read is a duplicate, False otherwise
    """
    return read["seq_id"] != read["prim_align_dup_exemplar"]


def count_reads_per_taxid(
    data: Iterator[dict[str, str]], taxid_field: str = "aligner_taxid_lca"
) -> tuple[Counter[TaxId], Counter[TaxId]]:
    """
    Count total and deduplicated reads per taxonomic ID.

    Args:
        data: Iterator of read records as dictionaries
        taxid_field: Field name containing the taxonomic ID

    Returns:
        Tuple of (total_counts, deduplicated_counts) as Counters
    """
    total = Counter()
    dedup = Counter()
    for read in data:
        taxid = read[taxid_field]
        total[taxid] += 1
        if not is_duplicate(read):
            dedup[taxid] += 1
    return total, dedup


def build_tree(
    tax_data: Iterator[dict[str, str]],
    child_field: str = "taxid",
    parent_field: str = "parent_taxid",
) -> Tree:
    """
    Build a taxonomic tree from taxonomy data.

    Args:
        tax_data: Iterator of taxonomy records as dictionaries
        child_field: Field name containing child taxonomic ID
        parent_field: Field name containing parent taxonomic ID

    Returns:
        Dictionary mapping parent IDs to lists of child IDs
    """
    tree = defaultdict(list)
    for taxon in tax_data:
        tree[taxon[parent_field]].append(taxon[child_field])
    return tree


def parents(tree: Tree) -> set[TaxId]:
    """Get all parent nodes from a tree."""
    return set(tree.keys())


def children(tree: Tree) -> set[TaxId]:
    """Get all child nodes from a tree."""
    return {child for lst in tree.values() for child in lst}


def nodes(tree: Tree) -> set[TaxId]:
    """Get all nodes (parents and children) from a tree."""
    return parents(tree) | children(tree)


def roots(tree: Tree) -> set[TaxId]:
    """Find root nodes of a tree (parent nodes that are not also child nodes)."""
    return parents(tree) - children(tree)


def aggregate_counts(node_counts: Counter[TaxId], tree: Tree) -> Counter[TaxId]:
    """
    Aggregate read counts for each clade (node and all its descendants).

    Args:
        node_counts: Counter of reads per taxonomic ID
        tree: Taxonomic tree structure

    Returns:
        Counter mapping taxonomic IDs to their clade counts
    """
    clade_counts = Counter()

    # Depth-first search of the tree, store results as you go
    def dfs(node: TaxId) -> int:
        # Start with this node's own count
        total = node_counts[node]

        # Add counts from all descendants
        for child in tree[node]:
            total += dfs(child)

        clade_counts[node] = total
        return total

    for root in roots(tree):
        dfs(root)

    return clade_counts


def write_output_tsv(
    output_path: str,
    tree: Tree,
    reads_total: Counter[TaxId],
    reads_dedup: Counter[TaxId],
    clade_counts_total: Counter[TaxId],
    clade_counts_dedup: Counter[TaxId],
) -> None:
    """
    Write taxonomic read counts to a TSV file.

    Args:
        output_path: Path to output TSV file
        tree: Taxonomic tree structure
        reads_total: Total read counts per taxonomic ID
        reads_dedup: Deduplicated read counts per taxonomic ID
        clade_counts_total: Total clade counts per taxonomic ID
        clade_counts_dedup: Deduplicated clade counts per taxonomic ID
    """
    with open_by_suffix(output_path, "w") as outfile:
        fieldnames = [
            "taxid",
            "parent_taxid",
            "reads_direct_total",
            "reads_direct_dedup",
            "reads_clade_total",
            "reads_clade_dedup",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        # Write rows in depth-first order
        def dfs(node: TaxId, parent = ".") -> None:
            row = {
                "taxid": node,
                "parent_taxid": parent,
                "reads_direct_total": reads_total[node],
                "reads_direct_dedup": reads_dedup[node],
                "reads_clade_total": clade_counts_total[node],
                "reads_clade_dedup": clade_counts_dedup[node],
            }
            # Only print clades that have some reads
            if row["reads_clade_total"] > 0:
                writer.writerow(row)
            for child in tree[node]:
                dfs(child, parent = node)

        for root in roots(tree):
            dfs(root)



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", help="Path to read TSV with LCA assignments.")
    parser.add_argument(
        "--taxdb", help="Path to taxonomy database with taxid and parent_taxid."
    )
    parser.add_argument("--output", help="Path to output TSV.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reads_total, reads_dedup = count_reads_per_taxid(read_tsv(args.reads))
    tree = build_tree(read_tsv(args.taxdb))
    clade_counts_total = aggregate_counts(reads_total, tree)
    clade_counts_dedup = aggregate_counts(reads_dedup, tree)
    write_output_tsv(
        args.output,
        tree,
        reads_total,
        reads_dedup,
        clade_counts_total,
        clade_counts_dedup,
    )


if __name__ == "__main__":
    main()
