#!/usr/bin/env python3

import argparse
import gzip
import json

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict


def process_fastq_file(filepath):
    """Process a single fastq file and return length distribution."""
    lengths = defaultdict(int)

    with gzip.open(filepath, "rt") as inf:
        for title, sequence, quality in FastqGeneralIterator(inf):
            seql = len(sequence)
            lengths[seql] += 1

    return dict(lengths)

def main():
    parser = argparse.ArgumentParser(
        description="Process fastq files and aggregate read lengths by barcode"
    )
    parser.add_argument("-i", dest="input_fastqc", help="Input fastqc file")
    parser.add_argument("-o", dest="output_json", help="Output JSON file path")
    parser.add_argument("-stage", dest="stage", help="Stage name")
    parser.add_argument("-sample", dest="sample", help="Sample name")
    args = parser.parse_args()

    filepath = args.input_fastqc
    # Process the file and get length distribution
    lengths = process_fastq_file(filepath)

    json_structure = {
        args.sample: {
            args.stage: lengths
        }
    }

    with open(args.output_json, "w") as outf:
        json.dump(json_structure, outf, indent=2, sort_keys=True)

if __name__ == "__main__":
    main()
