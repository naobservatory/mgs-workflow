#! /usr/bin/env python

import os
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--batch_size", type=int, default=1024**3)  # 1GiB
parser.add_argument("--pod5_dir", type=str)
parser.add_argument("--output_dir", type=str)
args = parser.parse_args()

batches = []
current_batch = []
current_batch_size = 0
batch_num = 0

for fname in os.listdir(args.pod5_dir):
    path = os.path.join(args.pod5_dir, fname)
    size = os.path.getsize(path)

    if current_batch_size + size > args.batch_size:
        batch_num += 1
        batch_dir = os.path.join(args.output_dir, f"batch_{batch_num:04d}")
        os.makedirs(batch_dir, exist_ok=True)
        for pod5_file in current_batch:
            shutil.copy(pod5_file, batch_dir)
        current_batch = []
        current_batch_size = 0

    current_batch.append(path)
    current_batch_size += size

# Handle the remaining files in the last batch
if current_batch:
    batch_num += 1
    batch_dir = os.path.join(args.output_dir, f"batch_{batch_num:04d}")
    os.makedirs(batch_dir, exist_ok=True)
    for pod5_file in current_batch:
        shutil.copy(pod5_file, batch_dir)