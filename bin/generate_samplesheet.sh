#!/bin/bash

##### Input parameters #####

# Initialize variables
dir_path=""
forward_suffix=""
reverse_suffix=""
s3=0
single_end=0
output_path="samplesheet.csv"  # Default output path

# Function to print usage
print_usage() {
    echo "Usage:"
    echo "For paired-end reads:"
    echo "  $0 --dir_path <path> --forward_suffix <suffix> --reverse_suffix <suffix> [--s3] [--output_path <path>]"
    echo "For single-end reads:"
    echo "  $0 --dir_path <path> --single_end [--s3] [--output_path <path>]"
    echo
    echo "Options:"
    echo "  --dir_path        Directory containing FASTQ files"
    echo "  --forward_suffix  Suffix for forward reads (required for paired-end only)"
    echo "  --reverse_suffix  Suffix for reverse reads (required for paired-end only)"
    echo "  --single_end     Flag for single-end data"
    echo "  --s3             Flag for S3 bucket access"
    echo "  --output_path    Output path for samplesheet (default: samplesheet.csv)"
}

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
        --single_end)
            single_end=1
            shift
            ;;
        --output_path)
            output_path="$2"
            shift 2
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Validate arguments based on single_end flag
if [ -z "$dir_path" ]; then
    echo "Error: dir_path is required."
    print_usage
    exit 1
fi

if [ $single_end -eq 0 ]; then
    # Paired-end validation
    if [[ -z "$forward_suffix" || -z "$reverse_suffix" ]]; then
        echo "Error: forward_suffix and reverse_suffix are required for paired-end reads."
        print_usage
        exit 1
    fi
fi

# Display the parameters
echo "Parameters:"
echo "dir_path: $dir_path"
echo "single_end: $single_end"
echo "s3: $s3"
echo "output_path: $output_path"
if [ $single_end -eq 0 ]; then
    echo "forward_suffix: $forward_suffix"
    echo "reverse_suffix: $reverse_suffix"
fi

##### Script #####

# Create header based on single_end flag
if [ $single_end -eq 0 ]; then
    echo "sample,fastq_1,fastq_2" > "$output_path"
else
    echo "sample,fastq" > "$output_path"
fi

# Ensure dir_path ends with a '/'
if [[ "$dir_path" != */ ]]; then
    dir_path="${dir_path}/"
fi

# Get file listing based on s3 flag
if [ $s3 -eq 1 ]; then
    listing=$(aws s3 ls ${dir_path} | awk '{print $4}')
else
    listing=$(ls ${dir_path} | awk '{print $1}')
fi

# Process files based on single_end flag
if [ $single_end -eq 0 ]; then
    # Paired-end processing
    echo "$listing" | grep "${forward_suffix}\.fastq\.gz$" | while read -r forward_read; do
        sample=$(echo "$forward_read" | sed -E "s/${forward_suffix}\.fastq\.gz$//")
        reverse_read=$(echo "$listing" | grep "${sample}${reverse_suffix}\.fastq\.gz$")
        if [ -n "$reverse_read" ]; then
            echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> "$output_path"
        fi
    done
else
    # Single-end processing - just process all fastq.gz files
    echo "$listing" | grep "\.fastq\.gz$" | while read -r read_file; do
        sample=$(echo "$read_file" | sed -E "s/\.fastq\.gz$//")
        echo "$sample,${dir_path}${read_file}" >> "$output_path"
    done
fi

echo "CSV file '$output_path' has been created."