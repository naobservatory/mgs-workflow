#!/usr/bin/env python

import os
import subprocess
import multiprocessing
import argparse
from pathlib import Path
import logging
import sys
import math
import glob
import re
import random

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def count_fasta_files(directory):
    return len(glob.glob(os.path.join(directory, "*.fna.gz")))

def download_chunk(chunk_info):
    chunk_file, output_dir, chunk_num, total_chunks = chunk_info
    with open(chunk_file, 'r') as f:
        taxid_count = sum(1 for _ in f)
    logging.info(f"Starting download for chunk {chunk_num+1}/{total_chunks} with {taxid_count} taxids")
    # Set up output directories and files
    chunk_output_dir = os.path.join(output_dir, f"genbank_genomes_chunk_{chunk_num}")
    metadata_file = f"chunk_{chunk_num}.metadata.txt"
    # Construct the ncbi-genome-download command
    cmd = [
        "ncbi-genome-download",
        "--section", "genbank",
        "--formats", "fasta",
        "--flat-output",
        "--taxids", chunk_file,
        "--metadata-table", metadata_file,
        "-o", chunk_output_dir,
        "viral"
    ]
    try:
        # Run the command and capture output
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        genome_count = count_fasta_files(chunk_output_dir)
        logging.info(f"Download complete for chunk {chunk_num+1}/{total_chunks}. Downloaded {genome_count} genomes.")
        return True, metadata_file, chunk_output_dir, None
    except subprocess.CalledProcessError as e:
        if "No downloads matched your filter" in e.stderr:
            error_msg = f"No matching downloads for chunk {chunk_num+1}/{total_chunks}"
            logging.warning(error_msg)
            return False, None, None, error_msg
        else:
            error_msg = f"Error processing chunk {chunk_num+1}/{total_chunks}: {e}\n"
            error_msg += f"STDOUT: {e.stdout}\n"
            error_msg += f"STDERR: {e.stderr}"
            logging.error(error_msg)
            return False, None, None, error_msg

def split_taxids(taxids_file, n_chunks):
    """Split taxid file into a specified number of chunks."""
    # Read all taxids from file
    with open(taxids_file, "r") as f:
        taxids = f.readlines()
    # Shuffle taxids to even workload across chunks
    random.shuffle(taxids)
    # Calculate chunk indices
    logging.info(f"Splitting {len(taxids)} taxids into {n_chunks} chunks.")
    chunk_size = math.ceil(len(taxids) / n_chunks)
    chunks = [taxids[i:i + chunk_size] for i in range(0, len(taxids), chunk_size)]
    chunk_files = []
    for i, chunk in enumerate(chunks, 1):
        chunk_file = f"taxid_chunk_{i}.txt"
        with open(chunk_file, 'w') as f:
            f.writelines(chunk)
        chunk_files.append(chunk_file)
    logging.info(f"Chunk splitting complete.")
    return(chunk_files)

def download_chunks(chunk_files, output_dir):
    """Download genomes from taxid chunks in parallel."""
    num_threads = len(chunk_files)
    logging.info(f"Starting genome download process with {num_threads} threads.")
    chunk_info = [(chunk_files[n], output_dir, n, num_threads) for n in range(num_threads)]
    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.map(download_chunk, chunk_info)
    logging.info("All downloads complete.")
    # Check for any errors
    errors = [error for _, _, _, error in results if error is not None]
    true_error = False
    if errors:
        for error in errors:
            if "No matching downloads" not in error:
                sys.exit(1)
        logging.warning("Some chunks had no matching downloads.")
    return results

def merge_chunks(results, output_dir, metadata_file):
    """Combine chunk output files into single output."""
    # Merge metadata
    logging.info("Merging metadata across chunks.")
    # Combine metadata from all successful downloads
    header = None
    data_rows = []
    for success, chunk_metadata_file, _, _ in results:
        if success and chunk_metadata_file:
            with open(chunk_metadata_file, 'r') as infile:
                chunk_lines = infile.readlines()
                if len(chunk_lines) > 0:
                    if header is None:
                        header = chunk_lines[0]
                    data_rows.extend(chunk_lines[1:])
    logging.info(f"{len(data_rows)} metadata rows found across chunks.")
    logging.info("Correcting metadata paths.")
    # Update paths in the chunk's data rows
    chunk_dir_pattern = re.escape(output_dir) + r"/genbank_genomes_chunk_\d+"
    #for i in range(len(data_rows)):
    #    data_rows[i] = re.sub(f"{chunk_dir_pattern}/", f"{output_dir}/", data_rows[i])
    with open(metadata_file, 'w') as outfile:
        if header:
            outfile.write(header)
            outfile.writelines(data_rows)
        else:
            logging.warning("No metadata found in any of the chunks.")
    logging.info("Metadata combined.")
    # Combine downloaded genomes into a single directory
    #for success, _, chunk_output_dir, _ in results:
    #    if success and chunk_output_dir:
    #        os.system(f"mv {chunk_output_dir}/* {output_dir}/")
    #logging.info("Genomes combined into single directory.")

def clean_chunks(chunk_files, results):
    """Clean up temporary chunk files."""
    logging.info("Cleaning up temporary files")
    for chunk_file in chunk_files:
        os.remove(chunk_file)
    for _, chunk_metadata_file, chunk_output_dir, _ in results:
        if chunk_metadata_file:
            os.remove(chunk_metadata_file)
        #if chunk_output_dir:
        #    os.rmdir(chunk_output_dir)
    logging.info("Cleaning complete.")

def main(taxids_file, num_threads, output_dir, metadata_file, random_seed):
    if random_seed is not None:
        random.seed(random_seed)
    # Generate chunk files
    chunk_files = split_taxids(taxids_file, num_threads)
    # Download chunk files
    os.makedirs(output_dir, exist_ok=True)
    chunk_results = download_chunks(chunk_files, output_dir)
    # Combine chunks
    merge_chunks(chunk_results, output_dir, metadata_file)
    # Clean temporary files
    clean_chunks(chunk_files, chunk_results)
    logging.info(f"Process complete. Results available in '{output_dir}' directory and '{metadata_file}'")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Download viral genomes in parallel")
    parser.add_argument("taxids_file", type=str, help="Path to the file containing taxids")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of parallel threads to use")
    parser.add_argument("--output", type=str, default="genbank_genomes",
                        help="Path to the output directory for downloaded genomes")
    parser.add_argument("--metadata", type=str, default="ncbi_fetch_metadata.txt",
                        help="Path to the output metadata file")
    parser.add_argument("--seed", type=int, help="Random seed for shuffling taxids")
    args = parser.parse_args()
    # Run the main function with provided arguments
    main(args.taxids_file, args.threads, args.output, args.metadata, args.seed)
