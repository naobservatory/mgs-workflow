#!/bin/bash


##### Input parameters #####

# Initialize variables
dir_path=""
forward_suffix=""
reverse_suffix=""
s3=0
output_path="samplesheet.csv"  # Default output path

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
            s3=1
            shift
            ;;
        --output_path)
            output_path="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if all required parameters are provided
if [[ -z "$dir_path" || -z "$forward_suffix" || -z "$reverse_suffix" ]]; then
    echo "Error: dir_path, forward_suffix, and reverse_suffix are required."
    echo "Usage: $0 --dir_path <path> --forward_suffix <suffix> --reverse_suffix <suffix> [--s3] [--output_path <path>]"
    exit 1
fi

# Display the parameters
echo "Parameters:"
echo "dir_path: $dir_path"
echo "forward_suffix: $forward_suffix"
echo "reverse_suffix: $reverse_suffix"
echo "s3: $s3"
echo "output_path: $output_path"

#### EXAMPLES ####

# dir_path="" # Cannot share this as it's restricted, but imagine the read looks like this 
# forward_suffix="_S[0-9]_L[0-9][0-9][0-9]_R1_001"
# reverse_suffix="_S[0-9]_L[0-9][0-9][0-9]_R2_001"
# s3=1

# Example 2
#  dir_path="s3://nao-mgs-public/PRJNA661613/raw/"
#  forward_suffix="_1"
#  reverse_suffix="_2"
#  s3=1

# Example 3
# dir_path="raw"
# forward_suffix="_1"
# reverse_suffix="_2"
# s3=0

##### Script #####

echo "sample,fastq_1,fastq_2" > "$output_path"

# Ensure dir_path ends with a '/'
if [[ "$dir_path" != */ ]]; then
    dir_path="${dir_path}/"
fi

listing=0

if [ $s3 -eq 1 ]; then
    listing=$(aws s3 ls ${dir_path} | awk '{print $4}')
else
    listing=$(ls ${dir_path} | awk '{print $1}')
fi


echo "$listing" | grep "${forward_suffix}\.fastq\.gz$" | while read -r forward_read; do
    sample=$(echo "$forward_read" | sed -E "s/${forward_suffix}\.fastq\.gz$//")
    reverse_read=$(echo "$listing" | grep "${sample}${reverse_suffix}\.fastq\.gz$")
    # If sample + reverse_suffix exists in s3_listing, then add to samplesheet
    if [ -n "$reverse_read" ]; then
        echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> "$output_path"
    fi
done

echo "CSV file '$output_path' has been created."
