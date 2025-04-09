#!/usr/bin/env python

#=======================================================================
# Preamble
#=======================================================================

# Import modules
import logging
import argparse
import pandas as pd
import time
from datetime import datetime, timezone

# Configure logging
class UTCFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime('%Y-%m-%d %H:%M:%S UTC')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter('[%(asctime)s] %(message)s')
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)

#=======================================================================
# Auxiliary functions
#=======================================================================

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process tabular output from VSEARCH clustering.")
    parser.add_argument("vsearch_db", help="Path to tabular output from VSEARCH clustering.")
    parser.add_argument("output_path", help="Output path for processed data frame.")
    args = parser.parse_args()
    logger.info("Command-line arguments parsed.")
    logger.info(f"VSEARCH DB path: {args.vsearch_db}")
    logger.info(f"Output path: {args.output_path}")
    return args.vsearch_db, args.output_path

def parse_vsearch_db(vsearch_path):
    """Parse VSEARCH clustering output into a DataFrame."""
    logger.info("Importing VSEARCH clustering output.")
    col_names = ["record_type", "cluster_id", "size", "percent_identity",
                 "orientation", "empty_0", "empty_1", "cigar",
                 "id", "cluster_rep_id"]
    vsearch_db = pd.read_csv(vsearch_path, sep="\t", header=None,
                             names=col_names, dtype=str)
    record_counts = vsearch_db["record_type"].value_counts()
    logger.info(
        f"Imported {record_counts.get('C', 0)} clusters, "
        f"{record_counts.get('S', 0)} representative sequences, "
        f" and {record_counts.get('H', 0)} non-representative hit sequences."
        )
    assert record_counts['C'] == record_counts['S'], \
            "Cluster and representative counts should match."
    return vsearch_db

def process_cluster_records(cluster_db):
    """Process cluster records into a structured format."""
    logger.info("Processing cluster records.")
    cluster_db = cluster_db.assign(
        cluster_size = cluster_db["size"].astype(int),
        cluster_rep_id = cluster_db["id"]
    )[["cluster_id", "cluster_size", "cluster_rep_id"]]
    logger.info(f"Processed {len(cluster_db)} cluster records.")
    logger.info(f"Output DB columns: {cluster_db.columns.tolist()}.")
    return cluster_db

def process_representative_records(rep_db):
    """Process cluster representative records into a structured format."""
    logger.info("Processing representative records.")
    rep_db = rep_db.assign(
        seq_length = rep_db["size"].astype(int),
        percent_identity = 100.0,
        orientation = "+",
        cigar = rep_db["size"].astype(str) + "M",
        seq_id = rep_db["id"],
        cluster_rep_id = rep_db["id"],
        is_cluster_rep = True
        )[["seq_id", "cluster_id", "cluster_rep_id", "seq_length",
            "is_cluster_rep", "percent_identity", "orientation", "cigar"]]
    logger.info(f"Processed {len(rep_db)} representative records.")
    logger.info(f"Output DB columns: {rep_db.columns.tolist()}.")
    return rep_db

def process_hit_records(hit_db):
    """Process hit records into a structured format."""
    logger.info("Processing hit records.")
    hit_db = hit_db.assign(
        seq_length = hit_db["size"].astype(int),
        percent_identity = hit_db["percent_identity"].astype(float),
        seq_id = hit_db["id"],
        is_cluster_rep = False
        )[["seq_id", "cluster_id", "cluster_rep_id", "seq_length",
            "is_cluster_rep", "percent_identity", "orientation", "cigar"]]
    logger.info(f"Processed {len(hit_db)} hit records.")
    logger.info(f"Output DB columns: {hit_db.columns.tolist()}.")
    return hit_db

def process_vsearch_db(vsearch_db):
    """Process VSEARCH tabular output into a structured format."""
    # First split DataFrame by record type
    logger.info("Processing VSEARCH DB records.")
    dbs_split = {label: db for label, db in vsearch_db.groupby("record_type")}
    assert "C" in dbs_split, "Cluster records should be present."
    assert "S" in dbs_split, "Representative records should be present."
    # Process cluster records
    cluster_db = process_cluster_records(dbs_split["C"].copy())
    rep_db = process_representative_records(dbs_split["S"].copy())
    if "H" in dbs_split:
        hit_db = process_hit_records(dbs_split["H"].copy())
    # Combine cluster and hit records
    if "H" in dbs_split:
        logger.info("Combining representative and hit records.")
        sequence_db = pd.concat([rep_db, hit_db], ignore_index=True)
    else:
        sequence_db = rep_db
    # Merge cluster information into sequence records
    logger.info("Merging cluster information into sequence records.")
    output_db = pd.merge(sequence_db, cluster_db,
                         on=["cluster_id", "cluster_rep_id"],
                         how="inner")
    assert len(output_db) == len(sequence_db), \
        "Output DB should match sequence DB length."
    # Output DB columns should equal sequence_db plus cluster columns
    assert output_db.shape[1] == sequence_db.shape[1] + cluster_db.shape[1] - 2, \
        "Output DB should have correct number of columns."
    logger.info(f"Done. Final output dimensions: {output_db.shape}.")
    logger.info(f"Output DB columns: {output_db.columns.tolist()}.")
    return output_db

def save_output(output_db, out_path):
    """Save processed DataFrame to output path."""
    logger.info(f"Saving output to {out_path}.")
    output_db.to_csv(out_path, sep="\t", index=False)
    logger.info("Output saved successfully.")

#=======================================================================
# Main function
#=======================================================================

def main():
    start_time = time.time()
    logger.info("Initializing process.")
    # Parse command-line arguments
    vsearch_path, out_path = parse_arguments()
    # Import VSEARCH clustering output (tab-delimited, no header)
    vsearch_db = parse_vsearch_db(vsearch_path)
    # Process VSEARCH DB into output format
    output_db = process_vsearch_db(vsearch_db)
    # Save output
    save_output(output_db, out_path)
    end_time = time.time()
    logger.info("Process completed.")
    logger.info("Total time elapsed: %.2f seconds" % (end_time - start_time))

if __name__ == "__main__":
    main()
