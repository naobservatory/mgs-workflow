#!/usr/bin/env python3

# Import modules
import sys
import argparse
import pandas as pd
import time
import datetime
import pysam
import gzip
import bz2
from collections import defaultdict
from Bio import SeqIO

# Utility functions

def print_log(message):
    print("[", datetime.datetime.now(), "]\t", message, sep="", file=sys.stderr)

def open_by_suffix(filename, mode="r"):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode)

def join_line(fields):
    "Convert a list of arguments into a TSV string for output."
    return("\t".join(map(str, fields)) + "\n")

# Alignment-level functions

def parse_sam_alignment(read, genbank_metadata, viral_taxids, clean_query_record):
    """Parse a Minimap2 SAM alignment."""
    out = {}
    out["query_name"] = read.query_name

    reference_genome_name = read.reference_name
    out["minimap2_genome_id_primary"] = reference_genome_name
    reference_taxid = extract_viral_taxid(reference_genome_name, genbank_metadata, viral_taxids)
    out["minimap2_taxid_primary"] = reference_taxid


    # Adding original read sequence and quality
    if read.is_reverse:
        # When minimap2 maps to the RC version of a strand, if returns the RC version of the read
        clean_query_record = clean_query_record.reverse_complement()

    query_seq_clean = clean_query_record.seq
    query_qual_clean_numeric = clean_query_record.letter_annotations["phred_quality"]
    # Turn numeric quality scores back into Phred+33 scores.
    query_qual_clean = "".join([chr(x+33) for x in query_qual_clean_numeric])
    query_len_clean = len(query_seq_clean)

    out["minimap2_read_length"] = read.query_length
    out["minimap2_map_qual"] = read.mapping_quality
    out["minimap2_ref_start"] = read.reference_start
    out["minimap2_ref_end"] = read.reference_end
    out["minimap2_alignment_start"] = read.query_alignment_start
    out["minimap2_alignment_end"] = read.query_alignment_end
    out["minimap2_cigar"] = read.cigarstring
    out["minimap2_edit_distance"] = read.get_tag("NM")
    out["minimap2_alignment_score"] = read.get_tag("AS")
    out["minimap2_query_sequence"] = read.query_sequence
    out["query_seq_clean"] = query_seq_clean
    out["query_qual_clean"] = query_qual_clean
    out["query_len_clean"] = query_len_clean

    return out


def extract_viral_taxid(genome_id, genbank_metadata, viral_taxids):
    """Extract taxid from the appropriate field of Genbank metadata."""
    try:
        taxid, species_taxid = genbank_metadata[genome_id]
        if taxid in viral_taxids:
            return taxid
        if species_taxid in viral_taxids:
            return species_taxid
        return taxid
    except KeyError:
        raise ValueError(f"No matching genome ID found: {genome_id}")


# File-level functions

def process_sam(sam_file, out_file, genbank_metadata, viral_taxids, clean_read_dict):
    """Process a Minimap2 SAM file."""
    with open_by_suffix(out_file, "w") as out_fh:
        header = (
            "query_name\t"
            "minimap2_genome_id_primary\t"
            "minimap2_taxid_primary\t"
            "minimap2_read_length\t"
            "minimap2_map_qual\t"
            "minimap2_ref_start\t"
            "minimap2_ref_end\t"
            "minimap2_alignment_start\t"
            "minimap2_alignment_end\t"
            "minimap2_cigar\t"
            "minimap2_edit_distance\t"
            "minimap2_alignment_score\t"
            "minimap2_query_sequence\t"
            "query_seq_clean\t"
            "query_qual_clean\t"
            "query_len_clean\n"
        )
        out_fh.write(header)
        try:
            with pysam.AlignmentFile(sam_file, "r") as sam_file:
                for read in sam_file:
                    print_log(f"Processing read: {read.query_name}")
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    read_id = read.query_name

                    clean_query_record = clean_read_dict[read_id]
                    line = parse_sam_alignment(read, genbank_metadata, viral_taxids, clean_query_record)
                    if line is None:
                        continue
                    line_keys = line.keys()
                    test_key_line = "\t".join(line_keys) + "\n"
                    assert test_key_line == header

                    out_fh.write(join_line(line.values()))

        except Exception as e:
            import traceback
            error_detail = traceback.format_exc()
            print_log(f"Error processing SAM file: {str(e)}")
            print_log(f"Error details: {error_detail}")
            raise

def parse_arguments():
    """Parse and return command-line arguments.

    Returns:
        Parsed argument namespace
    """
    parser = argparse.ArgumentParser(
        description="Process Minimap2 SAM output into a TSV with viral alignment information."
    )
    parser.add_argument(
        "-a", "--sam",
        type=str,
        required=True,
        help="Path to Minimap2 SAM alignment file."
    )
    parser.add_argument(
        "-r", "--reads",
        type=lambda f: open_by_suffix(f, "r"),
        default=sys.stdin,
        help="Path to FASTQ that contains the non-masked version of HV reads (default: stdin)."
    )
    parser.add_argument(
        "-m", "--metadata",
        required=True,
        help="Path to Genbank metadata file containing genomeID and taxid information."
    )
    parser.add_argument(
        "-v", "--viral_db",
        required=True,
        help="Path to TSV containing viral taxonomic information."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output path for processed data frame."
    )
    return parser.parse_args()

def main():
    # Parse arguments
    args = parse_arguments()

    try:
        sam_file = args.sam
        clean_reads = args.reads
        meta_path = args.metadata
        vdb_path = args.viral_db

        out_file = args.output
        # Start time tracking
        print_log("Starting process.")
        start_time = time.time()
        # Print parameters
        print_log("SAM file path: {}".format(sam_file))
        print_log("Genbank metadata file path: {}".format(meta_path))
        print_log("Viral DB file path: {}".format(vdb_path))
        print_log("Output path: {}".format(out_file))

        # Import metadata and viral DB
        print_log("Importing Genbank metadata file...")
        meta_db = pd.read_csv(meta_path, sep="\t", dtype=str)
        genbank_metadata = {
            genome_id: [taxid, species_taxid]
            for genome_id, taxid, species_taxid in zip(meta_db["genome_id"], meta_db["taxid"], meta_db["species_taxid"])
        }
        print_log("Importing viral DB file...")
        virus_db = pd.read_csv(vdb_path, sep="\t", dtype=str)
        viral_taxids = set(virus_db["taxid"].values)
        print_log(f"Virus DB imported. {len(virus_db)} total viral taxids.")


        print_log(f"Imported {len(viral_taxids)} virus taxa.")

        # Import clean reads
        clean_read_dict = {}
        for record in SeqIO.parse(clean_reads, "fastq"):
            read_id = record.id
            clean_read_dict[read_id] = record

        # Process SAM
        print_log("Processing SAM file...")
        process_sam(sam_file, out_file, genbank_metadata, viral_taxids, clean_read_dict)
        print_log("File processed.")

        # Finish time tracking
        end_time = time.time()
        print_log(f"Total time elapsed: {end_time - start_time:.2f} seconds")

    except Exception as e:
        print_log(f"Error: {str(e)}")
        sys.exit(1)
    finally:
        args.reads.close()  # Only close args.reads since args.sam is a string


if __name__ == "__main__":
    main()


