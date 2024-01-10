#!/usr/bin/env python

# Import modules
import re
import sys
import argparse
import pandas as pd
from collections import Counter

def get_parents(nodes_path):
    """Convert a nodes DMP file into a dictionary of child:parent taxon mappings."""
    parents = {}  # child_taxid -> parent_taxid
    with open(nodes_path) as inf:
        for line in inf:
            child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
            child_taxid = int(child_taxid)
            parent_taxid = int(parent_taxid)
            parents[child_taxid] = parent_taxid
    return(parents)

def get_assignments(kraken_path, parents):
    """Count direct and indirect clade assignments in a Kraken2 output file."""
    direct_assignments = Counter()
    clade_assignments = Counter()
    with open(kraken_path) as inf:
        for line in inf:
            line = line.strip()
            if not line: continue
            _, _, name_and_taxid, _, encoded_hits = line.split("\t")
            taxid, = re.findall("^.*[(]taxid ([0-9]+)[)]$", name_and_taxid)
            taxid = int(taxid)
            direct_assignments[taxid] += 1
            while True:
                clade_assignments[taxid] += 1
                if taxid in [0, 1]: break
                taxid = parents[taxid]
    df_clade = pd.DataFrame(clade_assignments, index=["clade_assignments"]).T.sort_index().reset_index()
    df_direct = pd.DataFrame(direct_assignments, index=["direct_assignments"]).T.sort_index().reset_index()
    df_out = df_clade.merge(df_direct, how = "outer").fillna(0).astype(int).rename(columns={"index":"taxid"})
    return(df_out)

def get_hits(kraken_path, parents):
    """Count direct and indirect partial hits in a Kraken2 output file."""
    direct_hits = Counter()
    clade_hits = Counter()
    with open(kraken_path) as inf:
        for line in inf:
            line = line.strip()
            if not line: continue
            _, _, name_and_taxid, _, encoded_hits = line.split("\t")
            direct_incremented = set()
            clade_incremented = set()
            for hit in re.findall("([0-9]+):", encoded_hits):
                hit = int(hit)
                if hit not in direct_incremented:
                    direct_hits[hit] += 1
                    direct_incremented.add(hit)
                while hit not in clade_incremented:
                    clade_hits[hit] += 1
                    clade_incremented.add(hit)
                    if hit in [0, 1]: break
                    hit = parents[hit]
    df_clade = pd.DataFrame(clade_hits, index=["clade_hits"]).T.sort_index().reset_index()
    df_direct = pd.DataFrame(direct_hits, index=["direct_hits"]).T.sort_index().reset_index()
    df_out = df_clade.merge(df_direct, how = "outer").fillna(0).astype(int).rename(columns={"index":"taxid"})
    return(df_out)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Count the number of reads matching each node in a taxonomic tree by Kraken2.")
    parser.add_argument("kraken_output", help="Path to Kraken2 output file.")
    parser.add_argument("nodes_dmp", help="Path to DMP file specifying tree structure.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    args = parser.parse_args()
    kraken_path = args.kraken_output
    nodes_path = args.nodes_dmp
    out_path = args.output_path
    # Extract tree structure
    parents = get_parents(nodes_path)
    # Get hits and assignments
    df_hits = get_hits(kraken_path, parents)
    df_assignments = get_assignments(kraken_path, parents)
    df_all = df_hits.merge(df_assignments, how = "outer").fillna(0).astype(int)
    # Write output
    gzip_output = True if out_path.split(".")[-1] == "gz" else False
    if gzip_output:
        df_all.to_csv(out_path, compression = "gzip", index = False)
    else:
        df_all.to_csv(out_path, index = False)

if __name__ == "__main__":
    main()

