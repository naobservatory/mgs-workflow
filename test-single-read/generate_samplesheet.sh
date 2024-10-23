#!/bin/bash


##### Input parameters #####

# Initialize variables
dir_path="s3://nao-mgs-simon/test_single_read/raw/"
forward_suffix=""
reverse_suffix=""
s3=1
single_end=1

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dir_path)
            dir_path="$2"
            shift 2
            ;;
        --forward_suffix)
            forward_suffix="$2"
            shift 2
            ;;
        --reverse_suffix)
            reverse_suffix="$2"
            shift 2
            ;;
        --s3)
            s3="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if all required parameters are provided
if [[ -z "$dir_path" || -z "$s3" ]]; then
    echo "Error: For both single-end and paired-end data, the directory path, single-end flag, and s3 flag are required."
    echo "Usage: $0 --dir_path <path> --single_end <0|1> --s3 <0|1>"
    exit 1
fi

# Check if suffixes are provided for paired-end data
if [ $single_end -eq 0 ] && [[ -z "$forward_suffix" || -z "$reverse_suffix" ]]; then
    echo "Error: Suffixes are required for paired-end data."
    echo "Usage: $0 --dir_path <path> --forward_suffix <suffix> --reverse_suffix <suffix> --s3 <0|1>"
    exit 1
fi

# Display the parameters
echo "Parameters:"
echo "dir_path: $dir_path"
echo "forward_suffix: $forward_suffix"
echo "reverse_suffix: $reverse_suffix"
echo "s3: $s3"
echo "single_end: $single_end"
#### EXAMPLES ####

# dir_path="" # Cannot share this as it's restricted, but imagine the read looks like this
# forward_suffix="_S[0-9]_L[0-9][0-9][0-9]_R1_001"
# reverse_suffix="_S[0-9]_L[0-9][0-9][0-9]_R2_001"
# s3=1
# single_end=0

# Example 2
#  dir_path="s3://nao-mgs-public/PRJNA661613/raw/"
#  forward_suffix="_1"
#  reverse_suffix="_2"
#  s3=1
# single_end=0

# Example 3
# dir_path="raw"
# forward_suffix="_1"
# reverse_suffix="_2"
# s3=0
# single_end=0

# Example 4
# dir_path="raw"
# forward_suffix=""
# reverse_suffix=""
# s3=0
# single_end=1

##### Script #####

## Single-end data

# Set the header based on single_end flag
if [ $single_end -eq 1 ]; then
    echo "sample,fastq" > samplesheet.csv
else
    echo "sample,fastq_1,fastq_2" > samplesheet.csv
fi

# Get file listing
if [ $s3 -eq 1 ]; then
    listing=$(aws s3 ls ${dir_path} | awk '{print $4}')
else
    listing=$(ls ${dir_path} | awk '{print $1}')
fi

echo "$listing" | grep "${forward_suffix}\.fastq\.gz$" | while read -r forward_read; do
    sample=$(echo "$forward_read" | sed -E "s/${forward_suffix}\.fastq\.gz$//")
    if [ $single_end -eq 1 ]; then
        echo "$sample,${dir_path}${forward_read}" >> samplesheet.csv
    else
        reverse_read=$(echo "$listing" | grep "${sample}${reverse_suffix}\.fastq\.gz$")
        # If sample + reverse_suffix exists in listing, then add to samplesheet
        if [ -n "$reverse_read" ]; then
            echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> samplesheet.csv
        fi
    fi
done

echo "CSV file 'samplesheet.csv' has been created."