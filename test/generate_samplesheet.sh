#!/bin/bash


##### Input parameters #####

# Path to the directory containing the fastq files
dir_path="raw"

# Suffix for forward reads of fastq files (e.g., "_1")
forward_suffix="_1"

# Suffix for reverse reads of fastq files (e.g., "_2")
reverse_suffix="_2"

# Middle part of the filename matched through regex; typically reads don't have this part (e.g., "_S1_L001" or "")
middle=""

# Set to 1 if the fastq files are in an S3 bucket, 0 if they are in a local directory
s3=0

##### Script #####

echo "sample,fastq_1,fastq_2" > samplesheet.csv

if [ $s3 -eq 1 ]; then
    s3_listing=$(aws s3 ls ${dir_path} | awk '{print $4}')
    
    echo "$s3_listing" | grep "${forward_suffix}\.fastq\.gz$" | while read -r forward_read; do
        reverse_read="${forward_read/${forward_suffix}/${reverse_suffix}}"
        if echo "$s3_listing" | grep -q "^${reverse_read}$"; then
            sample=$(echo "$forward_read" | sed -E "s/(.*)${forward_suffix}\.fastq\.gz/\1/" | sed -E "s/${middle}$//")
            echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> samplesheet.csv
        fi
    done
else
    find "$dir_path" -type f -name "*${forward_suffix}.fastq.gz" | while read -r forward_read; do
        reverse_read="${forward_read/${forward_suffix}/${reverse_suffix}}"
        if [ -f "$reverse_read" ]; then
            filename=$(basename "$forward_read")
            sample=$(echo "$filename" | sed -E "s/(.*)${forward_suffix}\.fastq\.gz/\1/" | sed -E "s/${middle}$//")
            echo "$sample,$forward_read,$reverse_read" >> samplesheet.csv
        fi
    done
fi

echo "CSV file 'samplesheet.csv' has been created."



#### EXAMPLES ####

# dir_path="" # Cannot share this as it's restricted, but imagine the read looks like this 
# forward_suffix="_R1_001"
# reverse_suffix="_R2_001"
# middle="_S[0-9]_L[0-9][0-9][0-9]"
# s3=1

# Example 2
#  dir_path="s3://nao-mgs-public/PRJNA661613/raw/"
#  forward_suffix="_1"
#  reverse_suffix="_2"
#  middle=""
#  s3=1

# Example 3
# dir_path="raw"
# forward_suffix="_1"
# reverse_suffix="_2"
# middle=""
# s3=0