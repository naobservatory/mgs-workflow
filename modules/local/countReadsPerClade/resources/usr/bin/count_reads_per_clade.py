import gzip
import csv
from collections import Counter, defaultdict
from sys import argv

def read_tsv(file_path, encoding='utf-8'):
    """
    Read a gzipped TSV file and yield rows one at a time.

    Args:
        file_path (str): Path to the gzipped TSV file
        encoding (str): File encoding (default: 'utf-8')

    Yields:
        dict: Dictionary representing each row
    """
    with gzip.open(file_path, 'rt', encoding=encoding) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            yield row

def is_duplicate(read):
    return read["seq_id"] == read["prim_align_dup_exemplar"]

def count_reads_per_taxid(data):
    total = Counter()
    dedup = Counter()
    for read in data:
        taxid = read["aligner_taxid_lca"]
        total[taxid] += 1
        if not is_duplicate(read):
            dedup[taxid] += 1
    return total, dedup

def build_tree(tax_data):
    # {parent: [child]}
    tree = defaultdict(list)
    for taxon in tax_data:
        tree[taxon["parent_taxid"]].append(taxon["taxid"])
    return tree

def roots(tree):
    # Find roots of a tree (parents nodes are not also child nodes)
    parents = set(tree.keys())
    children = {child for lst in tree.values() for child in lst} 
    return parents - children

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


def main():
    read_table = argv[1]
    tax_db = argv[2]
    reads_total, reads_dedup = count_reads_per_taxid(read_tsv(read_table))
    tree = build_tree(read_tsv(tax_db))
    clade_counts_total = aggregate_counts(reads_total, tree)
    clade_counts_dedup = aggregate_counts(reads_dedup, tree)

if __name__ == "__main__":
    main()
