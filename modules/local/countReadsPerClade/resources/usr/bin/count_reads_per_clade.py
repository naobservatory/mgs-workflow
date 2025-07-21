import argparse
import gzip
import bz2
import csv
from collections import Counter, defaultdict


def open_by_suffix(filename, mode="r", debug=False):
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode + "t")
    elif filename.endswith(".bz2"):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


def read_tsv(file_path):
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


def is_duplicate(read):
    return read["seq_id"] != read["prim_align_dup_exemplar"]


def count_reads_per_taxid(data, taxid_field="aligner_taxid_lca"):
    total = Counter()
    dedup = Counter()
    for read in data:
        taxid = read[taxid_field]
        total[taxid] += 1
        if not is_duplicate(read):
            dedup[taxid] += 1
    return total, dedup


def build_tree(tax_data, child_field="taxid", parent_field="parent_taxid"):
    # {parent: [child]}
    tree = defaultdict(list)
    for taxon in tax_data:
        tree[taxon[parent_field]].append(taxon[child_field])
    return tree


def parents(tree):
    return set(tree.keys())


def children(tree):
    return {child for lst in tree.values() for child in lst}


def nodes(tree):
    return parents(tree) | children(tree)


def roots(tree):
    # Find roots of a tree (parents nodes that are not also child nodes)
    return parents(tree) - children(tree)


def aggregate_counts(node_counts, tree):
    clade_counts = {}

    # Depth-first search of the tree, store results as you go
    def dfs(node):
        # Start with this node's own count
        total = node_counts.get(node, 0)

        # Add counts from all descendants
        for child in tree[node]:
            total += dfs(child)

        clade_counts[node] = total
        return total

    for root in roots(tree):
        dfs(root)

    return clade_counts


def write_output_tsv(
    output_path,
    tree,
    reads_total,
    reads_dedup,
    clade_counts_total,
    clade_counts_dedup,
):
    with open_by_suffix(output_path, "w") as outfile:
        fieldnames = [
            "taxid",
            "reads_direct_total",
            "reads_direct_dedup",
            "reads_clade_total",
            "reads_clade_dedup",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for node in nodes(tree):
            row = {
                "taxid": node,
                "reads_direct_total": reads_total.get(node, 0),
                "reads_direct_dedup": reads_dedup.get(node, 0),
                "reads_clade_total": clade_counts_total.get(node, 0),
                "reads_clade_dedup": clade_counts_dedup.get(node, 0),
            }
            # Only print clades that have some reads
            if row["reads_clade_total"] > 0:
                writer.writerow(row)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", help="Path to read TSV with LCA assignments.")
    parser.add_argument(
        "--taxdb", help="Path to taxonomy database with taxid and parent_taxid."
    )
    parser.add_argument("--output", help="Path to output TSV.")
    return parser.parse_args()


def main():
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
