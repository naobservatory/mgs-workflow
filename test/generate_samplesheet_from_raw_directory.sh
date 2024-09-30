#!/bin/bash

[[ $# -eq 0 || ($# -eq 1 && $1 == "-s") ]] && { echo "Usage: $0 [-s] <directory>"; exit 1; }

s3=false; [[ $1 == "-s" ]] && { s3=true; shift; }
dir=$1; out="samplesheet.csv"

echo "sample,fastq_1,fastq_2" > "$out"

if $s3; then
    aws s3 ls "$dir" | awk '{print $4}' | grep '_1\.fastq\.gz$' | while read -r file; do
        prefix=${file%_1.fastq.gz}
        echo "$prefix,$dir$file,$dir${prefix}_2.fastq.gz"
    done >> "$out"
else
    ls "$dir" | grep '_1\.fastq\.gz$' | while read -r file; do
        prefix=${file%_1.fastq.gz}
        echo "$prefix,$dir/$file,$dir/${prefix}_2.fastq.gz"
    done >> "$out"
fi

echo "Generated: $out"
