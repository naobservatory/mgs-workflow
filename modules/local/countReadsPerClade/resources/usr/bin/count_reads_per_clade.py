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
    # {parent: [children]}
    children = defaultdict(list)
    for taxon in tax_data:
        children[taxon["parent_taxid"]].append(taxon["taxid"])
    return children

def get_roots(children):
    # Find roots of a tree
    all_children = {child for lst in children.values() for child in lst} 
    parents = set(children.keys())
    return parents - all_children

def aggregate_counts(node_counts, children):
    roots = get_roots(children)

    # Post-order traversal to compute clade counts
    clade_counts = {}

    def dfs(node):
        # Start with this node's own count
        total = node_counts.get(node, 0)

        # Add counts from all descendants
        for child in children[node]:
            total += dfs(child)

        clade_counts[node] = total
        return total

    for root in roots:
        dfs(root)

    return clade_counts

def main():
    read_table = read_tsv(argv[1])
    tax_db = read_tsv(argv[2])
    reads_total, reads_dedup = count_reads_per_taxid(read_table)
    children = build_tree(tax_db)
    clade_counts = aggregate_counts(reads_total, children)
    # Note: total should equal the sum of the roots
    for taxid, count in clade_counts.items():
        if count > 0:
            print(taxid, count, sep="\t")
    print(clade_counts["1"])
    print(reads_total.total())

if __name__ == "__main__":
    main()
