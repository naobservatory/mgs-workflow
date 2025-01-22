#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Generate a samplesheet CSV from FASTQ files in a directory or S3 bucket."
    )

    # Required / optional arguments
    parser.add_argument("--dir_path", required=True, help="Directory or S3 path containing FASTQ files")
    parser.add_argument("--single_end", action="store_true",
                        help="Flag indicating data is single-end (omit for paired-end)")
    parser.add_argument("--forward_suffix", default="",
                        help="Regex for forward reads, e.g. '_R1_001' (required if not single-end)")
    parser.add_argument("--reverse_suffix", default="",
                        help="Regex for reverse reads, e.g. '_R2_001' (required if not single-end)")
    parser.add_argument("--s3", action="store_true",
                        help="Flag indicating files are on an S3 bucket")
    parser.add_argument("--output_path", default="samplesheet.csv",
                        help="Output CSV path (default: samplesheet.csv)")
    parser.add_argument("--group_file", default="",
                        help="Optional CSV with columns 'sample,group' to merge grouping info")
    parser.add_argument("--group_across_illumina_lanes", action="store_true",
                        help="Group samples across Illumina lanes by removing '_Lnnn' from sample name")
    
    args = parser.parse_args()

    # Validate arguments
    if not args.single_end:
        # Paired-end requires forward_suffix and reverse_suffix
        if not args.forward_suffix or not args.reverse_suffix:
            parser.error("--forward_suffix and --reverse_suffix are required for paired-end data.")

    if args.single_end:
        if args.forward_suffix or args.reverse_suffix:
            parser.error("--forward_suffix and --reverse_suffix are not allowed for single-end data.")

    if "s3://" in args.dir_path and not args.s3:
        parser.error("--s3 is required when dir_path is an S3 bucket.")

    if args.group_file and args.group_across_illumina_lanes:
        parser.error("Provide at most one of --group_file and --group_across_illumina_lanes.")

    # Print parameters (optional, for clarity)
    print("Parameters:")
    print(f"  dir_path: {args.dir_path}")
    print(f"  single_end: {args.single_end}")
    print(f"  s3: {args.s3}")
    print(f"  output_path: {args.output_path}")
    print(f"  group_file: {args.group_file}")
    print(f"  group_across_illumina_lanes: {args.group_across_illumina_lanes}")
    if not args.single_end:
        print(f"  forward_suffix: {args.forward_suffix}")
        print(f"  reverse_suffix: {args.reverse_suffix}")
    print()

    # Ensure dir_path ends with '/'
    if not args.dir_path.endswith("/"):
        args.dir_path += "/"

    # 1. Get the file listing
    listing = get_file_listing(args.dir_path, args.s3)

    # 2. Build a dictionary of samples -> list of fastq paths
    #    For paired-end, each sample should have 2 files: forward and reverse.
    #    For single-end, each sample should have 1 file.
    samples_dict = build_samples_dict(
        listing=listing,
        dir_path=args.dir_path,
        single_end=args.single_end,
        forward_suffix=args.forward_suffix,
        reverse_suffix=args.reverse_suffix
    )

    # 3. Incorporate grouping info
    #    Optionally from a group_file or by trimming `_Lnnn`.
    groups_dict = {}
    if args.group_file:
        # Read group file into a dict
        groups_dict = read_group_file(args.group_file)
    elif args.group_across_illumina_lanes:
        # We'll generate group info by removing _Lnnn from sample name
        for sample_name in samples_dict:
            groups_dict[sample_name] = re.sub(r"_L\d{3}$", "", sample_name)

    # 4. Write to CSV
    write_samplesheet(
        samples_dict=samples_dict,
        groups_dict=groups_dict,
        output_path=args.output_path,
        single_end=args.single_end
    )

    print(f"\nDone. CSV file '{args.output_path}' has been created.")


def get_file_listing(dir_path, is_s3):
    """
    Returns a list of filenames from either an S3 bucket or a local directory.
    """
    if is_s3:
        # Use aws s3 ls
        try:
            cmd = ["aws", "s3", "ls", dir_path]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True, text=True)
            # Output lines look like: "YYYY-MM-DD HH:MM:SS size filename"
            # We parse out the last token (filename)
            listing = [line.split()[-1] for line in result.stdout.strip().split("\n") if line.strip()]
            return listing
        except subprocess.CalledProcessError as e:
            print(f"Error listing S3 bucket: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Local directory
        if not os.path.isdir(dir_path):
            print(f"Error: {dir_path} is not a directory or cannot be accessed.", file=sys.stderr)
            sys.exit(1)
        return os.listdir(dir_path)


def build_samples_dict(listing, dir_path, single_end, forward_suffix, reverse_suffix):
    """
    Builds and returns a dictionary:
      sample_name -> (for single-end)  [fq_file]
                     (for paired-end)  [fq_file_fwd, fq_file_rev]
    """
    # Pre-compile regex for performance if suffix is used
    # We expect the files to end with ".fastq.gz" (the script assumes that)
    fastq_pattern = re.compile(r"\.fastq\.gz$")
    
    if not single_end:
        # Paired-end
        fwd_regex = re.compile(forward_suffix + r"\.fastq\.gz$")
        rev_regex = re.compile(reverse_suffix + r"\.fastq\.gz$")

        # We'll maintain a dict: sample_name -> {"fwd": None, "rev": None}
        samples_dict = {}

        for filename in listing:
            # Ignore anything not .fastq.gz
            if not fastq_pattern.search(filename):
                continue

            full_path = dir_path + filename

            # If it matches the forward suffix
            if fwd_regex.search(filename):
                # Derive sample name by removing the forward suffix
                sample_name = fwd_regex.sub("", filename)
                if sample_name not in samples_dict:
                    samples_dict[sample_name] = {"fwd": None, "rev": None}
                samples_dict[sample_name]["fwd"] = full_path

            elif rev_regex.search(filename):
                sample_name = rev_regex.sub("", filename)
                if sample_name not in samples_dict:
                    samples_dict[sample_name] = {"fwd": None, "rev": None}
                samples_dict[sample_name]["rev"] = full_path
            else:
                # It's a .fastq.gz but doesn't match either suffix—ignore
                continue

        # Convert that structure to final dict of sample_name -> [fwd, rev]
        # Excluding samples that don't have both reads
        final_dict = {}
        for sample, fr_dict in samples_dict.items():
            if fr_dict["fwd"] and fr_dict["rev"]:
                final_dict[sample] = [fr_dict["fwd"], fr_dict["rev"]]
        return final_dict

    else:
        # Single-end
        samples_dict = {}
        for filename in listing:
            if fastq_pattern.search(filename):
                sample_name = fastq_pattern.sub("", filename)
                samples_dict[sample_name] = [dir_path + filename]
        return samples_dict


def read_group_file(group_file_path):
    """
    Reads a CSV group file with at least two columns: sample,group
    Returns a dict: {sample_name -> group_name}
    """
    groups_dict = {}
    with open(group_file_path, "r", newline="") as gf:
        reader = csv.reader(gf)
        header = next(reader, None)
        if not header:
            print("Group file is empty or missing a header.", file=sys.stderr)
            sys.exit(1)
        # We expect the first two columns to be sample,group
        for row in reader:
            if len(row) >= 2:
                sample_name = row[0]
                group_name = row[1]
                groups_dict[sample_name] = group_name
    return groups_dict


def write_samplesheet(samples_dict, groups_dict, output_path, single_end):
    """
    Writes out a CSV file with either:
      single_end -> sample,fastq[,group]
      paired_end -> sample,fastq_1,fastq_2[,group]
    """
    # Decide on header
    if single_end:
        header = ["sample", "fastq"]
    else:
        header = ["sample", "fastq_1", "fastq_2"]

    # If we have at least one group, append 'group' to header
    # We determine if grouping is in effect by checking whether `groups_dict` is non-empty.
    # (If you specifically want 'NA' for everything if no group, you can do that too.)
    has_groups = len(groups_dict) > 0
    if has_groups:
        header.append("group")

    # Write rows
    with open(output_path, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)

        for sample_name, fq_files in samples_dict.items():
            # For single-end, fq_files has length 1
            # For paired-end, fq_files has length 2
            if single_end:
                row = [sample_name, fq_files[0]]
            else:
                row = [sample_name, fq_files[0], fq_files[1]]

            if has_groups:
                group_val = groups_dict.get(sample_name, "NA")
                row.append(group_val)

            writer.writerow(row)


if __name__ == "__main__":
    main()