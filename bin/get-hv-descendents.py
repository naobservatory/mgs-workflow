#!/usr/bin/env python

import argparse
from collections import defaultdict
import subprocess
import json
import multiprocessing
import sys
import time

def run_gimme_taxa(taxid):
    # Run gimme_taxa.py for a given taxid
    cmd = ["gimme_taxa.py", str(taxid)]
    try:
        p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return p.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching descendants for taxid {taxid}: {e}")
        print(f"Error output: {e.stderr}")
        raise  # Re-raise the exception if it's not about missing database

def fetch_descendents(hv_taxid):
    print(f"\tFetching descendants of {hv_taxid}", flush=True)
    descendents = [hv_taxid]
    try:
        output = run_gimme_taxa(hv_taxid)
    except Exception as e:
        print(f"Error: {e}")
        raise  # Re-raise the exception to be handled by the calling function
    desc_split = output.split("\n")
    # Add descendent taxids to lists
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
    print(f"- {len(descendents) - 1} total descendents fetched for {hv_taxid}.")
    return hv_taxid, descendents

def main():
    parser = argparse.ArgumentParser(description="Get descendents of human virus taxids")
    parser.add_argument("input_file", help="Path to the input file containing human virus taxids")
    parser.add_argument("output_list", help="Path to the output file for all taxids")
    parser.add_argument("output_json", help="Path to the output JSON file for taxid descendents")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), 
                        help="Number of threads to use (default: number of CPU cores)")
    args = parser.parse_args()
    # Load human viruses
    print("Importing virus db...", end="")
    human_viruses = {}
    with open(args.input_file) as inf:
        for line in inf:
            taxid, name = line.strip().split("\t")
            human_viruses[int(taxid)] = name
    print("done.")
    # Boot NCBI database through an unthreaded trial run
    print("Booting NCBI database...", end="")
    cmd = ["gimme_taxa.py", "662620"]
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print("done.")
    # Get descendents of each viral taxid
    print("Fetching viral descendents:")
    with multiprocessing.Pool(processes=args.threads) as pool:
        try:
            results = pool.map(fetch_descendents, human_viruses.keys())
        except Exception as e:
            print(f"Error during descendents fetching: {e}")
            sys.exit(1)
    taxid_descendents = dict(results)
    taxid_unique = set(taxid for descendents in taxid_descendents.values() for taxid in descendents)
    print("Fetching complete.")
    # Write taxid list to text file
    with open(args.output_list, "w") as outf:
        for taxid in sorted(taxid_unique):
            outf.write(f"{taxid}\n")
    # Write ancestor:descendent mappings to JSON
    with open(args.output_json, "w") as outf:
        json.dump(taxid_descendents, outf)

if __name__ == "__main__":
    main()
