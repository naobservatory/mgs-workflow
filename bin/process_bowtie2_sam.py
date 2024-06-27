#!/usr/bin/env python

# Import modules
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
    print("[", datetime.datetime.now(), "]\t", message, sep="")

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
        print(query_value)
        print(fields)
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

def process_sam_alignment(sam_line, genomeid_taxid_map, paired):
    """Process a SAM alignment line."""
    fields_in = sam_line.strip().split("\t")
    out = {}
    out["query_name"] = fields_in[0]
    out.update(process_sam_flags(fields_in[1]))
    out["genome_id"] = fields_in[2]
    out["taxid"] = int(genomeid_taxid_map[fields_in[2]][0])
    # TODO: Get taxid from genome_id
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

def get_line(query_name, genome_id, taxid, fragment_length, best_alignment_score_fwd, best_alignment_score_ref,
             next_alignment_score_fwd, next_alignment_score_rev, edit_distance_fwd, edit_distance_rev,
             ref_start_fwd, ref_start_rev, map_qual_fwd, map_qual_rev, cigar_fwd, cigar_rev, query_len_fwd, query_len_rev,
             query_seq_fwd, query_seq_rev, pair_status):
    """Convert a list of arguments into an output line."""
    fields = [query_name, genome_id, taxid, fragment_length, best_alignment_score_fwd, best_alignment_score_ref,
              next_alignment_score_fwd, next_alignment_score_rev, edit_distance_fwd, edit_distance_rev,
              ref_start_fwd, ref_start_rev, map_qual_fwd, map_qual_rev, cigar_fwd, cigar_rev, query_len_fwd, query_len_rev,
              query_seq_fwd, query_seq_rev, pair_status]
    fields_joined = join_line(fields)
    return(fields_joined)

def process_sam_unpaired_pair(read_dict):
    """Process an unpaired alignment from a paired SAM file."""
    # Specify unpaired attributes
    if read_dict["is_mate_1"]:
        best_alignment_score_fwd, best_alignment_score_rev = read_dict["alignment_score"], "NA"
        next_alignment_score_fwd, next_alignment_score_rev = read_dict["next_best_alignment"], "NA"
        edit_distance_fwd, edit_distance_rev = read_dict["edit_distance"], "NA"
        ref_start_fwd, ref_start_rev = read_dict["ref_start"], "NA"
        map_qual_fwd, map_qual_rev = read_dict["map_qual"], "NA"
        cigar_fwd, cigar_rev = read_dict["cigar"], "NA"
        query_len_fwd, query_len_rev = read_dict["query_len"], "NA"
        query_seq_fwd, query_seq_rev = read_dict["query_seq"], "NA"
    else:
        best_alignment_score_fwd, best_alignment_score_rev = "NA", read_dict["alignment_score"]
        next_alignment_score_fwd, next_alignment_score_rev = "NA", read_dict["next_best_alignment"]
        edit_distance_fwd, edit_distance_rev = "NA", read_dict["edit_distance"]
        ref_start_fwd, ref_start_rev = "NA", read_dict["ref_start"]
        map_qual_fwd, map_qual_rev = "NA", read_dict["map_qual"]
        cigar_fwd, cigar_rev = "NA", read_dict["cigar"]
        query_len_fwd, query_len_rev = "NA", read_dict["query_len"]
        query_seq_fwd, query_seq_rev = "NA", read_dict["query_seq"]
    return get_line(read_dict["query_name"], read_dict["genome_id"], read_dict["taxid"], read_dict["fragment_length"],
                    best_alignment_score_fwd, best_alignment_score_rev, next_alignment_score_fwd, next_alignment_score_rev,
                    edit_distance_fwd, edit_distance_rev, ref_start_fwd, ref_start_rev, map_qual_fwd, map_qual_rev,
                    cigar_fwd, cigar_rev, query_len_fwd, query_len_rev, query_seq_fwd, query_seq_rev, read_dict["pair_status"])

def line_from_valid_pair(fwd_dict, rev_dict):
    """Generate an output line from a validated (concordant or discordant) pair of SAM alignments."""
    query_name = fwd_dict["query_name"]
    genome_id = fwd_dict["genome_id"]
    taxid = fwd_dict["taxid"]
    fragment_length = fwd_dict["fragment_length"]
    best_alignment_score_fwd = fwd_dict["alignment_score"]
    best_alignment_score_rev = rev_dict["alignment_score"]
    next_alignment_score_fwd = fwd_dict["next_best_alignment"]
    next_alignment_score_rev = rev_dict["next_best_alignment"]
    edit_distance_fwd = fwd_dict["edit_distance"]
    edit_distance_rev = rev_dict["edit_distance"]
    ref_start_fwd = fwd_dict["ref_start"]
    ref_start_rev = rev_dict["ref_start"]
    map_qual_fwd = fwd_dict["map_qual"]
    map_qual_rev = rev_dict["map_qual"]
    cigar_fwd = fwd_dict["cigar"]
    cigar_rev = rev_dict["cigar"]
    query_len_fwd = fwd_dict["query_len"]
    query_len_rev = rev_dict["query_len"]
    query_seq_fwd = fwd_dict["query_seq"]
    query_seq_rev = rev_dict["query_seq"]
    pair_status = fwd_dict["pair_status"]
    return get_line(query_name, genome_id, taxid, fragment_length,
                    best_alignment_score_fwd, best_alignment_score_rev, next_alignment_score_fwd, next_alignment_score_rev,
                    edit_distance_fwd, edit_distance_rev, ref_start_fwd, ref_start_rev, map_qual_fwd, map_qual_rev,
                    cigar_fwd, cigar_rev, query_len_fwd, query_len_rev, query_seq_fwd, query_seq_rev, pair_status)

def process_sam_concordant_pair(fwd_dict, rev_dict):
    """Process a concordant pair of SAM alignments."""
    # Verify concordant alignment
    try:
        assert fwd_dict["query_name"] == rev_dict["query_name"]
        assert fwd_dict["genome_id"] == rev_dict["genome_id"]
        assert fwd_dict["pair_status"] == rev_dict["pair_status"] == "CP"
        assert fwd_dict["fragment_length"] == rev_dict["fragment_length"]
        assert fwd_dict["ref_start"] == rev_dict["mate_ref_start"]
        assert rev_dict["ref_start"] == fwd_dict["mate_ref_start"]
        assert fwd_dict["in_pair"] == rev_dict["in_pair"] == True
        assert fwd_dict["proper_paired_alignment"] == rev_dict["proper_paired_alignment"] == True
    except AssertionError:
        print(fwd_line)
        print(rev_line)
        raise
    # Generate output line
    return line_from_valid_pair(fwd_dict, rev_dict)

def process_sam_discordant_pair(fwd_dict, rev_dict):
    """Process a discordant pair of SAM alignments."""
    # Check if both reads align to the same genome.
    if fwd_dict["genome_id"] != rev_dict["genome_id"]: # If not, process as separate alignments
        return [process_sam_unpaired_pair(fwd_dict), process_sam_unpaired_pair(rev_dict)]
    # Verify discordant same-subject alignment
    try:
        assert fwd_dict["query_name"] == rev_dict["query_name"]
        assert fwd_dict["genome_id"] == rev_dict["genome_id"]
        assert fwd_dict["pair_status"] == rev_dict["pair_status"] == "DP"
        assert fwd_dict["fragment_length"] == rev_dict["fragment_length"]
        assert fwd_dict["ref_start"] == rev_dict["mate_ref_start"]
        assert rev_dict["ref_start"] == fwd_dict["mate_ref_start"]
        assert fwd_dict["in_pair"] == rev_dict["in_pair"] == True
    except AssertionError:
        print(fwd_line)
        print(rev_line)
        raise
    return [line_from_valid_pair(fwd_dict, rev_dict)]

# File-level functions
def write_sam_headers_paired(out_file):
    """Write header line to new TSV."""
    headers = ["query_name", "genome_id", "taxid", "fragment_length", "best_alignment_score_fwd", "best_alignment_score_rev", "next_alignment_score_fwd", "next_alignment_score_rev", "edit_distance_fwd", "edit_distance_rev", "ref_start_fwd", "ref_start_rev", "map_qual_fwd", "map_qual_rev", "cigar_fwd", "cigar_rev", "query_len_fwd", "query_len_rev", "query_seq_fwd", "query_seq_rev", "pair_status"]
    header_line = join_line(headers)
    out_file.write(header_line)
    return None

def process_sam_alignments_paired(fwd_line, rev_line, in_file, out_file, genomeid_taxid_map):
    """Run through a paired SAM file & process into TSV."""
    pr=functools.partial(process_sam_alignments_paired, in_file=in_file, out_file=out_file, genomeid_taxid_map=genomeid_taxid_map)
    if fwd_line is not None: 
        fwd_dict = process_sam_alignment(fwd_line, genomeid_taxid_map, True)
        if fwd_dict["pair_status"] not in ["UU", "UP", "CP", "DP"]:
            raise ValueError("Invalid pair status: {}, {}".format(fwd_dict["query_name"], fwd_dict["pair_status"]))
        elif fwd_dict["pair_status"] == "UU":
            raise ValueError("Unpaired read in paired SAM file: {}".format(fwd_dict["query_name"]))
        elif fwd_dict["pair_status"] == "UP":
            print_log("Unpaired pair: {}".format(fwd_dict["query_name"]))
            line = process_sam_unpaired_pair(fwd_dict) # Process forward read, then move on one line and retry
            out_file.write(line)
            return(pr(rev_line, get_next_alignment(in_file)))
        else:
            if rev_line is None:
                raise ValueError("Forward read is paired but reverse read is missing: {}".format(fwd_dict["query_name"]))
            rev_dict = process_sam_alignment(rev_line, genomeid_taxid_map, True)
            if rev_dict["query_name"] != fwd_dict["query_name"]:
                msg = "Alignments should be paired but read IDs mismatch: {}, {}"
                raise ValueError(msg.format(fwd_dict["query_name"], rev_dict["query_name"]))
            elif fwd_dict["pair_status"] == "CP": # Process reads as a pair, then move on two lines
                line = process_sam_concordant_pair(fwd_dict, rev_dict)
                out_file.write(line)
            elif fwd_dict["pair_status"] == "DP":
                lines = process_sam_discordant_pair(fwd_dict, rev_dict)
                for line in lines:
                    out_file.write(line)
            return(pr(get_next_alignment(in_file), get_next_alignment(in_file)))
    elif rev_line is not None:
        rev_dict = process_sam_alignment(rev_line, genomeid_taxid_map, True)
        msg = "Invalid data: reverse line exists without forward line: {}".format(rev_dict["query_name"])
        raise ValueError(msg)
    else:
        return None

def process_paired_sam(sam_path, out_path, genomeid_taxid_map):
    """Process paired SAM file into a TSV."""
    with open_by_suffix(sam_path) as inf, open_by_suffix(out_path, "w") as outf:
        head = write_sam_headers_paired(outf) # Write headers
        # Process file & generate content
        fwd_line = get_next_alignment(inf)
        rev_line = get_next_alignment(inf)
        while True:
            if fwd_line is not None: # If forward line exists
                fwd_dict = process_sam_alignment(fwd_line, genomeid_taxid_map, True)
                if fwd_dict["pair_status"] not in ["UU", "UP", "CP", "DP"]:
                    raise ValueError("Invalid pair status: {}, {}".format(fwd_dict["query_name"], fwd_dict["pair_status"]))
                elif fwd_dict["pair_status"] == "UU":
                    raise ValueError("Unpaired read in paired SAM file: {}".format(fwd_dict["query_name"]))
                elif fwd_dict["pair_status"] == "UP":
                    print_log("Unpaired pair: {}".format(fwd_dict["query_name"]))
                    line = process_sam_unpaired_pair(fwd_dict) # Process forward read, then move on one line and retry
                    outf.write(line)
                    fwd_line = rev_line
                    rev_line = get_next_alignment(inf)
                    continue
                else:
                    if rev_line is None:
                        raise ValueError("Forward read is paired but reverse read is missing: {}".format(fwd_dict["query_name"]))
                    rev_dict = process_sam_alignment(rev_line, genomeid_taxid_map, True)
                    if rev_dict["query_name"] != fwd_dict["query_name"]:
                        msg = "Alignments should be paired but read IDs mismatch: {}, {}"
                        raise ValueError(msg.format(fwd_dict["query_name"], rev_dict["query_name"]))
                    elif fwd_dict["pair_status"] == "CP": # Process reads as a pair, then move on two lines
                        line = process_sam_concordant_pair(fwd_dict, rev_dict)
                        outf.write(line)
                    elif fwd_dict["pair_status"] == "DP":
                        lines = process_sam_discordant_pair(fwd_dict, rev_dict)
                        for line in lines:
                            outf.write(line)
                    fwd_line = get_next_alignment(inf)
                    rev_line = get_next_alignment(inf)
                    continue
            elif rev_line is not None:
                rev_dict = process_sam_alignment(rev_line, genomeid_taxid_map, True)
                msg = "Invalid data: reverse line exists without forward line: {}".format(rev_dict["query_name"])
                raise ValueError(msg)
            else:
                break

# TODO: Process unpaired SAM

# Main function
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process Bowtie2 SAM output into a TSV with additional information.")
    parser.add_argument("sam", help="Path to Bowtie2 SAM file.")
    parser.add_argument("genomeid_taxid_map", help="Path to JSON file containing genomeID/taxID mapping.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    parser.set_defaults(paired=True)
    parser.add_argument("-p", "--paired", dest="paired", action="store_true", help="Processed SAM file as containing paired read alignments (default: True).")
    parser.add_argument("-u", "--unpaired", dest="paired", action="store_false", help="Process SAM file as containing unpaired read alignments (default: False).")
    args = parser.parse_args()
    sam_path = args.sam
    out_path = args.output_path
    mapping_path = args.genomeid_taxid_map
    paired = args.paired
    # Start time tracking
    print_log("Starting process.")
    start_time = time.time()
    # Print parameters
    print_log("SAM file path: {}".format(sam_path))
    print_log("Mapping file path: {}".format(mapping_path))
    print_log("Output path: {}".format(out_path))
    print_log("Processing file as paired: {}".format(paired))
    # Import genomeID/taxID mapping
    print_log("Importing JSON mapping file...")
    with open_by_suffix(mapping_path) as inf:
        mapping = json.load(inf)
    print_log("JSON imported.")
    # Process SAM
    print_log("Processing SAM file...")
    if paired:
        process_paired_sam(sam_path, out_path, mapping)
    else:
        process_unpaired_sam(sam_path, out_path, mapping)
    print_log("File processed.")
    # Finish time tracking
    end_time = time.time()
    print_log("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
