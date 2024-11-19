#!/usr/bin/env python

# Import modules
import argparse
import time
import datetime
import gzip
import bz2
from Bio import SeqIO
from Bio import Seq
import os

def print_log(message):
    print("[", datetime.datetime.now(), "]\t", message, sep="")

def open_by_suffix(filename, mode="r"):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

def join_paired_reads(file1, file2, output_file, gap="N"):
    """Join non-overlapping paired-end reads from two FASTQ files."""
    with open_by_suffix(file1) as f1, open_by_suffix(file2) as f2, open_by_suffix(output_file, "w") as out:
        r1 = SeqIO.parse(f1, "fastq")
        r2 = SeqIO.parse(f2, "fastq")
        r  = zip(r1,r2)
        for fwd,rev in r:
            s1, q1 = str(fwd.seq), fwd.letter_annotations["phred_quality"]
            s2, q2 = str(rev.seq.reverse_complement()), rev.letter_annotations["phred_quality"][::-1]
            joined_seq = Seq.Seq(s1 + gap + s2)
            # Modify the description to include "joined" without duplicating the ID
            description_parts = fwd.description.split(maxsplit=1)
            if len(description_parts) > 1 and description_parts[0] == fwd.id:
                new_description = f"joined {description_parts[1]}"
            else:
                new_description = f"joined {fwd.description}"
            joined_qual = {"phred_quality": q1 + [0]*len(gap) + q2}
            joined_read = SeqIO.SeqRecord(joined_seq, id=fwd.id, name=fwd.name, 
                                          description=new_description,
                                          letter_annotations=joined_qual)
            SeqIO.write(joined_read, out, "fastq")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Join non-overlapping read pairs into single sequences.")
    parser.add_argument("fwd_reads", help="Path to FASTQ file containing forward reads.")
    parser.add_argument("rev_reads", help="Path to FASTQ file containing reverse reads.")
    parser.add_argument("output_file", help="Path to output FASTQ file containing joined reads.")
    parser.add_argument("-g", "--gap", help="Gap sequence separating joined reads. (Default: N)",
                        default="N", nargs=1)
    args = parser.parse_args()
    fwd_path = args.fwd_reads
    rev_path = args.rev_reads
    out_path = args.output_file
    gap = args.gap
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("Forward reads file: {}".format(fwd_path))
    print_log("Reverse reads file: {}".format(rev_path))
    print_log("Output reads file: {}".format(out_path))
    print_log("Joining string: {}".format(gap))
    # Run joining function
    print_log("Joining reads...")
    join_paired_reads(fwd_path, rev_path, out_path, gap) # NB: Currently implemented serially. Will investigate paralellising if necessary.
    print_log("...done.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
