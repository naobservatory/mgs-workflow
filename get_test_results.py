#! /usr/bin/env python3

import subprocess
import os
mgs_results_dir = "/home/ec2-user/code/mgs-workflow/mgs-results"

single_read_results_dir = "s3://nao-mgs-simon/test_single_read/output/"
single_read_work_dir = "s3://nao-mgs-simon/test_single_read/work/"
paired_read_results_dir = "s3://nao-mgs-simon/test_paired_end/output/"
paired_read_work_dir = "s3://nao-mgs-simon/test_paired_end/work/"

def get_test_results():
    mgs_results_single = os.path.join(mgs_results_dir, "single_end", "output")
    mgs_results_paired = os.path.join(mgs_results_dir, "paired_end", "output")

    os.makedirs(mgs_results_single, exist_ok=True)
    os.makedirs(mgs_results_paired, exist_ok=True)

    subprocess.run(["aws", "s3", "sync", single_read_results_dir, mgs_results_single])
    subprocess.run(["aws", "s3", "sync", paired_read_results_dir, mgs_results_paired])

def get_test_work_dirs():
    mgs_work_single = os.path.join(mgs_results_dir, "single_end", "work")
    mgs_work_paired = os.path.join(mgs_results_dir, "paired_end", "work")

    os.makedirs(mgs_work_single, exist_ok=True)
    os.makedirs(mgs_work_paired, exist_ok=True)

    subprocess.run(["aws", "s3", "sync", single_read_work_dir, mgs_work_single])
    subprocess.run(["aws", "s3", "sync", paired_read_work_dir, mgs_work_paired])

if __name__ == "__main__":
    get_test_results()
    get_test_work_dirs()