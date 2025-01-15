#!/usr/bin/env python

# Import modules
import sys
import re
import argparse
import pandas as pd
import time
import datetime
import gzip
import bz2
import json
import functools

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

# Single-line functions

def get_next_alignment(sam_file):
    """Iterate through SAM file lines until gets an alignment line, then returns."""
    while True:
        l = next(sam_file, "EOF") # Get next line
        if not l: continue # Skip empty lines
        if l.startswith("@"): continue # Skip header lines
        if l == "EOF": return None
        return(l)

def check_flag(flag_sum, flag_dict, flag_name, flag_value):
    """Check if a flag sum includes a specific flag and return adjusted flag sum."""
    flag_sum = int(flag_sum)
    if flag_sum >= flag_value:
        flag_dict[flag_name] = True
        flag_sum -= flag_value
    else:
        flag_dict[flag_name] = False
    return flag_sum, flag_dict

def process_sam_flags(flag_sum):
    """Extract individual flags from flag sum."""
    flag_dict = {}
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "is_mate_2", 128)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "is_mate_1", 64)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "mate_aligned_reverse", 32)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "aligned_reverse", 16)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "no_paired_alignments", 8)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "no_single_alignments", 4)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "proper_paired_alignment", 2)
    flag_sum, flag_dict = check_flag(flag_sum, flag_dict, "in_pair", 1)
    return(flag_dict)

def extract_option(opt_list, query_value, default=None):
    """Extract specific optional field from list of optional fields."""
    fields = [f for f in opt_list if query_value in f]
    if len(fields) == 0: return(default)
    try:
        assert len(fields) == 1
    except AssertionError:
        print_log(query_value)
        print_log(fields)
        raise
    field = fields[0]
    pattern = "{}:.:(.*)".format(query_value)
    out_value = re.findall(pattern, field)[0]
    return(out_value)

def extract_optional_fields(opt_list):
    """Extract relevant optional fields from Bowtie2 SAM alignment."""
    out = {}
    out["alignment_score"] = extract_option(opt_list, "AS")
    out["next_best_alignment"] = extract_option(opt_list, "XS")
    out["edit_distance"] = extract_option(opt_list, "NM")
    out["pair_status"] = extract_option(opt_list, "YT")
    out["mate_alignment_score"] = extract_option(opt_list, "YS")
    return(out)

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
        msg = "No matching genome ID found: {}".format(genome_id)
        raise ValueError(msg)

def process_sam_alignment(sam_line, genbank_metadata, viral_taxids, paired):
    """Process a SAM alignment line."""
    fields_in = sam_line.strip().split("\t")
    out = {}
    out["query_name"] = fields_in[0]
    out.update(process_sam_flags(fields_in[1]))
    out["genome_id"] = fields_in[2]
    out["taxid"] = int(extract_viral_taxid(fields_in[2], genbank_metadata, viral_taxids))
    out["ref_start"] = int(fields_in[3]) - 1 # Convert from 1-indexing to 0-indexing
    out["map_qual"] = int(fields_in[4])
    out["cigar"] = fields_in[5]
    if paired:
        out["mate_genome_id"] = fields_in[6]
        out["mate_ref_start"] = int(fields_in[7]) - 1 # Convert as above
        out["fragment_length"] = abs(int(fields_in[8]))
    out["query_seq"] = fields_in[9]
    out["query_len"] = len(fields_in[9])
    out["query_qual"] = fields_in[10]
    out.update(extract_optional_fields(fields_in[11:]))
    return(out)

def join_line(fields):
    "Convert a list of arguments into a TSV string for output."
    return("\t".join(map(str, fields)) + "\n")

def get_line(fields_dict):
    """Convert a dictionary of arguments into an output line."""
    fields = [fields_dict["query_name"],
              fields_dict["genome_id_fwd"], fields_dict["genome_id_rev"],
              fields_dict["taxid_fwd"], fields_dict["taxid_rev"],
              fields_dict["fragment_length_fwd"], fields_dict["fragment_length_rev"],
              fields_dict["best_alignment_score_fwd"], fields_dict["best_alignment_score_rev"],
              fields_dict["next_alignment_score_fwd"], fields_dict["next_alignment_score_rev"],
              fields_dict["edit_distance_fwd"], fields_dict["edit_distance_rev"],
              fields_dict["ref_start_fwd"], fields_dict["ref_start_rev"],
              fields_dict["map_qual_fwd"], fields_dict["map_qual_rev"],
              fields_dict["cigar_fwd"], fields_dict["cigar_rev"],
              fields_dict["query_len_fwd"], fields_dict["query_len_rev"],
              fields_dict["query_seq_fwd"], fields_dict["query_seq_rev"],
              fields_dict["query_qual_fwd"], fields_dict["query_qual_rev"],
              fields_dict["pair_status"]]
    fields_joined = join_line(fields)
    return(fields_joined)

def line_from_single(read_dict):
    """Generate an output line from a single SAM alignment dictionary."""
    out_dict = {"query_name": read_dict["query_name"], "pair_status": read_dict["pair_status"]}
    if read_dict["is_mate_1"]:
        out_dict["genome_id_fwd"], out_dict["genome_id_rev"] = read_dict["genome_id"], "NA"
        out_dict["taxid_fwd"], out_dict["taxid_rev"] = read_dict["taxid"], "NA"
        out_dict["fragment_length_fwd"], out_dict["fragment_length_rev"] = read_dict["fragment_length"], "NA"
        out_dict["best_alignment_score_fwd"], out_dict["best_alignment_score_rev"] = read_dict["alignment_score"], "NA"
        out_dict["next_alignment_score_fwd"], out_dict["next_alignment_score_rev"] = read_dict["next_best_alignment"], "NA"
        out_dict["edit_distance_fwd"], out_dict["edit_distance_rev"] = read_dict["edit_distance"], "NA"
        out_dict["ref_start_fwd"], out_dict["ref_start_rev"] = read_dict["ref_start"], "NA"
        out_dict["map_qual_fwd"], out_dict["map_qual_rev"] = read_dict["map_qual"], "NA"
        out_dict["cigar_fwd"], out_dict["cigar_rev"] = read_dict["cigar"], "NA"
        out_dict["query_len_fwd"], out_dict["query_len_rev"] = read_dict["query_len"], "NA"
        out_dict["query_seq_fwd"], out_dict["query_seq_rev"] = read_dict["query_seq"], "NA"
        out_dict["query_qual_fwd"], out_dict["query_qual_rev"] = read_dict["query_qual"], "NA"
    else:
        out_dict["genome_id_fwd"], out_dict["genome_id_rev"] = "NA", read_dict["genome_id"]
        out_dict["taxid_fwd"], out_dict["taxid_rev"] = "NA", read_dict["taxid"]
        out_dict["fragment_length_fwd"], out_dict["fragment_length_rev"] = "NA", read_dict["fragment_length"]
        out_dict["best_alignment_score_fwd"], out_dict["best_alignment_score_rev"] = "NA", read_dict["alignment_score"]
        out_dict["next_alignment_score_fwd"], out_dict["next_alignment_score_rev"] = "NA", read_dict["next_best_alignment"]
        out_dict["edit_distance_fwd"], out_dict["edit_distance_rev"] = "NA", read_dict["edit_distance"]
        out_dict["ref_start_fwd"], out_dict["ref_start_rev"] = "NA", read_dict["ref_start"]
        out_dict["map_qual_fwd"], out_dict["map_qual_rev"] = "NA", read_dict["map_qual"]
        out_dict["cigar_fwd"], out_dict["cigar_rev"] = "NA", read_dict["cigar"]
        out_dict["query_len_fwd"], out_dict["query_len_rev"] = "NA", read_dict["query_len"]
        out_dict["query_seq_fwd"], out_dict["query_seq_rev"] = "NA", read_dict["query_seq"]
        out_dict["query_qual_fwd"], out_dict["query_qual_rev"] = "NA", read_dict["query_qual"]
    return get_line(out_dict)

def line_from_pair(dict_1, dict_2):
    """Generate an output line from two SAM alignment dictionaries."""
    # Identify forward and reverse reads
    if dict_1["is_mate_1"] and dict_2["is_mate_1"]:
        raise ValueError("Both reads are forward reads: {}".format(dict_1["query_name"]))
    if not dict_1["is_mate_1"] and not dict_2["is_mate_1"]:
        raise ValueError("Both reads are reverse reads: {}".format(dict_1["query_name"]))
    fwd_dict = dict_1 if dict_1["is_mate_1"] else dict_2
    rev_dict = dict_1 if not dict_1["is_mate_1"] else dict_2
    # Prepare dictionary for output
    out_dict = {
        "query_name": fwd_dict["query_name"],
        "genome_id_fwd": fwd_dict["genome_id"],
        "genome_id_rev": rev_dict["genome_id"],
        "taxid_fwd": fwd_dict["taxid"],
        "taxid_rev": rev_dict["taxid"],
        "fragment_length_fwd": fwd_dict["fragment_length"],
        "fragment_length_rev": rev_dict["fragment_length"],
        "best_alignment_score_fwd": fwd_dict["alignment_score"],
        "best_alignment_score_rev": rev_dict["alignment_score"],
        "next_alignment_score_fwd": fwd_dict["next_best_alignment"],
        "next_alignment_score_rev": rev_dict["next_best_alignment"],
        "edit_distance_fwd": fwd_dict["edit_distance"],
        "edit_distance_rev": rev_dict["edit_distance"],
        "ref_start_fwd": fwd_dict["ref_start"],
        "ref_start_rev": rev_dict["ref_start"],
        "map_qual_fwd": fwd_dict["map_qual"],
        "map_qual_rev": rev_dict["map_qual"],
        "cigar_fwd": fwd_dict["cigar"],
        "cigar_rev": rev_dict["cigar"],
        "query_len_fwd": fwd_dict["query_len"],
        "query_len_rev": rev_dict["query_len"],
        "query_seq_fwd": fwd_dict["query_seq"],
        "query_seq_rev": rev_dict["query_seq"],
        "query_qual_fwd": fwd_dict["query_qual"],
        "query_qual_rev": rev_dict["query_qual"],
        "pair_status": fwd_dict["pair_status"]
        }
    return get_line(out_dict)

def check_pair_status(line_dict, paired):
    """Check if pair status is valid for SAM entry line."""
    pair_status = line_dict["pair_status"]
    if pair_status not in ["UU", "UP", "CP", "DP"]:
        raise ValueError(f"Invalid pair status: {line_dict['query_name']}, {pair_status}")
    elif paired and pair_status == "UU":
        raise ValueError(f"Unpaired read in paired SAM file: {line_dict['query_name']}")
    elif not paired and pair_status != "UU":
        raise ValueError(f"Paired read in unpaired SAM file: {line_dict['query_name']}")

# File-level functions
def write_sam_headers_paired(out_file):
    """Write header line to new TSV."""
    headers = ["query_name",
               "genome_id_fwd", "genome_id_rev",
               "taxid_fwd", "taxid_rev",
               "fragment_length_fwd", "fragment_length_rev",
               "best_alignment_score_fwd", "best_alignment_score_rev",
               "next_alignment_score_fwd", "next_alignment_score_rev",
               "edit_distance_fwd", "edit_distance_rev",
               "ref_start_fwd", "ref_start_rev",
               "map_qual_fwd", "map_qual_rev",
               "cigar_fwd", "cigar_rev",
               "query_len_fwd", "query_len_rev",
               "query_seq_fwd", "query_seq_rev",
               "query_qual_fwd", "query_qual_rev",
               "pair_status"]
    header_line = join_line(headers)
    out_file.write(header_line)
    return None

def process_paired_sam(inf, outf, genbank_metadata, viral_taxids):
    """Process paired SAM file into a TSV."""
    # Write headers
    head = write_sam_headers_paired(outf)
    # Process file and generate content
    fwd_line = get_next_alignment(inf)
    rev_line = get_next_alignment(inf)
    while True:
        if fwd_line is None: # If no line, check for reverse line then break
            if rev_line is not None: # Break if reverse line exists without forward line
                rev_dict = process_sam_alignment(rev_line, genbank_metadata, viral_taxids, True)
                msg = f"Invalid data: reverse line exists without forward line: {rev_dict['query_name']}"
                raise ValueError(msg)
            break
        # Extract forward read information and check pair status
        fwd_dict = process_sam_alignment(fwd_line, genbank_metadata, viral_taxids, True)
        check_pair_status(fwd_dict, True)
        if rev_line is None: # If no reverse line, check pair status, then process forward line as unpaired
            if fwd_dict["pair_status"] != "UP":
                msg = f"Forward read is paired but reverse read is missing: {fwd_dict['query_name']}"
                raise ValueError(msg)
            line = line_from_single(fwd_dict)
            outf.write(line)
            fwd_line = rev_line
            rev_line = get_next_alignment(inf)
            continue
        # Extract reverse read information and check pair status
        rev_dict = process_sam_alignment(rev_line, genbank_metadata, viral_taxids, True)
        check_pair_status(rev_dict, True)
        # Check for sorting
        if fwd_dict["query_name"] > rev_dict["query_name"]:
            msg = f"Reads are not sorted: encountered {fwd_dict['query_name']} before {rev_dict['query_name']}"
            raise ValueError(msg)
        # Check if read IDs match
        if rev_dict["query_name"] != fwd_dict["query_name"]:
            # If IDs mismatch, forward read should be unpaired
            if fwd_dict["pair_status"] != "UP":
                msg = f"Forward read is paired but reverse read is missing: {fwd_dict['query_name']}"
                raise ValueError(msg)
            line = line_from_single(fwd_dict)
            outf.write(line)
            fwd_line = rev_line
            rev_line = get_next_alignment(inf)
            continue
        # Check that pair statuses match and are valid
        if rev_dict["pair_status"] != fwd_dict["pair_status"]:
            msg = f"Pair status mismatch: {fwd_dict['query_name']}, {fwd_dict['pair_status']}, {rev_dict['pair_status']}"
            raise ValueError(msg)
        if fwd_dict["pair_status"] == "UP" or rev_dict["pair_status"] == "UP":
            msg = f"Both mates align but alignment is unpaired: {fwd_dict['query_name']}"
            raise ValueError(msg)
        # Process pair together
        line = line_from_pair(fwd_dict, rev_dict)
        outf.write(line)
        fwd_line = get_next_alignment(inf)
        rev_line = get_next_alignment(inf)
        continue

def parse_arguments():
    """Parse and return command-line arguments."""
    parser = argparse.ArgumentParser(description="Process a sorted Bowtie2 SAM output into a TSV with additional information.")
    parser.add_argument("-s", "--sam", type=lambda f: open_by_suffix(f, "r"), default=sys.stdin,
                        help="Path to Bowtie2 SAM file (default: stdin).")
    parser.add_argument("-m", "--metadata", required=True,
                        help="Path to Genbank download metadata file containing genomeID and taxid information.")
    parser.add_argument("-v", "--viral_db", required=True,
                        help="Path to TSV containing viral taxonomic information.")
    parser.add_argument("-o", "--output", type=lambda f: open_by_suffix(f, "w"), default=sys.stdout,
                        help="Output path for processed data frame (default: stdout).")
    parser.add_argument("-p", "--paired", dest="paired", action="store_true", help="Processed SAM file as containing paired read alignments (default: True).")
    parser.add_argument("-u", "--unpaired", dest="paired", action="store_false", help="Process SAM file as containing unpaired read alignments (default: False).")
    parser.set_defaults(paired=True)
    args = parser.parse_args()
    return args

# Main function
def main():
    # Parse arguments
    args = parse_arguments()
    try:
        sam_file = args.sam
        out_file = args.output
        meta_path = args.metadata
        vdb_path = args.viral_db
        paired = args.paired
        # Start time tracking
        print_log("Starting process.")
        start_time = time.time()
        # Print parameters
        print_log("SAM file path: {}".format(sam_file))
        print_log("Genbank metadata file path: {}".format(meta_path))
        print_log("Viral DB file path: {}".format(vdb_path))
        print_log("Output path: {}".format(out_file))
        print_log("Processing file as paired: {}".format(paired))
        # Import metadata and viral DB
        print_log("Importing Genmank metadata file...")
        meta_db = pd.read_csv(meta_path, sep="\t", dtype=str)
        gid_taxid_dict = {genome_id: [taxid, species_taxid]
                          for genome_id, taxid, species_taxid in
                          zip(meta_db["genome_id"],meta_db["taxid"],
                              meta_db["species_taxid"])}
        print_log("Importing viral DB file...")
        virus_db = pd.read_csv(vdb_path, sep="\t", dtype=str)
        virus_taxa = set(virus_db["taxid"].values)
        print_log("Imported {} virus taxa.".format(len(virus_taxa)))
        # Process SAM
        print_log("Processing SAM file...")
        if paired:
            process_paired_sam(sam_file, out_file, gid_taxid_dict, virus_taxa)
        else:
            process_unpaired_sam(sam_file, out_file, gid_taxid_dict, virus_taxa)
        print_log("File processed.")
        # Finish time tracking
        end_time = time.time()
        print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))
    finally:
        args.sam.close()
        args.output.close()

if __name__ == "__main__":
    main()
