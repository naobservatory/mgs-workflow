#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Merge JSON files")
    parser.add_argument("-i", dest="input_jsons", help="Input JSON files", nargs="+")
    parser.add_argument("-o", dest="output_json", help="Output JSON file")
    args = parser.parse_args()

    input_jsons = args.input_jsons
    output_json = args.output_json

    # Create output JSON structure
    lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for file in input_jsons:
        with open(file) as f:
            data = json.load(f)
            for sample, stages in data.items():
                for stage, items in stages.items():
                    for length, count in items.items():
                        lengths[sample][stage][length] = lengths[sample][stage][length] + count

    with open(output_json, "w") as f:
        json.dump(lengths, f, indent=2, sort_keys=True)

if __name__ == "__main__":
    main()
