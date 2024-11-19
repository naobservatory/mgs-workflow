#!/bin/bash

##### Input parameters #####

# Initialize variables
dir_path=""
forward_suffix=""
reverse_suffix=""
s3=0
single_end=0
output_path="samplesheet.csv"  # Default output path
group_file=""  # Optional parameter for the group file


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
        --group_file)  # Optional group file
            group_file="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Check if all required parameters are provided
if [[ -z "$dir_path" || -z "$single_end" ]]; then
    echo "Error: dir_path and single_end are required."
    echo -e "\nUsage: $0 [options]"
    echo -e "\nRequired arguments:"
    echo -e "  --dir_path <path>         Directory containing FASTQ files"
    echo -e "  --single_end              Flag for single-end data"
    echo -e "\nOptional arguments:"
    echo -e "  --forward_suffix <suffix>  When single_end is 0, suffix identifying forward reads, supports regex (e.g., '_R1_001' or '_1')"
    echo -e "  --reverse_suffix <suffix>  When single_end is 0, suffix identifying reverse reads, supports regex (e.g., '_R2_001' or '_2')"
    echo -e "  --s3                      Use if files are stored in S3 bucket"
    echo -e "  --output_path <path>      Output path for samplesheet [default: samplesheet.csv]"
    echo -e "  --group_file <path>       Path to group file for sample grouping [header column must have the names 'sample,group' in that order; additional columns may be included, however they will be ignored by the script]"
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
echo "group_file: $group_file"
if [ $single_end -eq 0 ]; then
    echo "forward_suffix: $forward_suffix"
    echo "reverse_suffix: $reverse_suffix"
fi


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

# Create a temporary file for the initial samplesheet
temp_samplesheet=$(mktemp)

# Create header based on single_end flag
if [ $single_end -eq 0 ]; then
    echo "sample,fastq_1,fastq_2" > "$temp_samplesheet"
else
    echo "sample,fastq" > "$temp_samplesheet"
fi
echo "group_file: $group_file"


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
        # If sample + reverse_suffix exists in s3_listing, then add to samplesheet
        if [ -n "$reverse_read" ]; then
            echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> "$temp_samplesheet"
        fi
    done
else
    # Single-end processing - just process all fastq.gz files
    echo "$listing" | grep "\.fastq\.gz$" | while read -r read_file; do
        sample=$(echo "$read_file" | sed -E "s/\.fastq\.gz$//")
        echo "$sample,${dir_path}${read_file}" >> "$temp_samplesheet"
    done
fi

# Check if group file is provided
if [[ -n "$group_file" ]]; then
    # Perform left join with group file
    awk -F',' 'NR==FNR{a[$1]=$2; next} FNR==1{print $0",group"} FNR>1{print $0","(a[$1]?a[$1]:"NA")}' "$group_file" "$temp_samplesheet" > "$output_path"
    echo "CSV file '$output_path' has been created with group information."
else
    # If no group file, just use the temporary samplesheet as the final output
    mv "$temp_samplesheet" "$output_path"
    echo "CSV file '$output_path' has been created without group information."
fi

# Remove temporary file if it still exists
[ -f "$temp_samplesheet" ] && rm "$temp_samplesheet"
