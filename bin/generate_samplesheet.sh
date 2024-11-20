#!/bin/bash

set -u
set -e

##### Input parameters #####

# Initialize variables
dir_path=""
forward_suffix=""
reverse_suffix=""
s3=0
output_path="samplesheet.csv"  # Default output path
group_file=""  # Optional parameter for the group file
group_across_illumina_lanes=false

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
        --group_file)  # Optional group file
            group_file="$2"
            shift 2
            ;;
        --group_across_illumina_lanes)
            group_across_illumina_lanes=true
            shift
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
    echo -e "\nUsage: $0 [options]"
    echo -e "\nRequired arguments:"
    echo -e "  --dir_path <path>         Directory containing FASTQ files"
    echo -e "  --forward_suffix <suffix>  Suffix identifying forward reads, supports regex (e.g., '_R1_001' or '_1')"
    echo -e "  --reverse_suffix <suffix>  Suffix identifying reverse reads, supports regex (e.g., '_R2_001' or '_2')"
    echo -e "\nOptional arguments:"
    echo -e "  --s3                      Use if files are stored in S3 bucket"
    echo -e "  --output_path <path>      Output path for samplesheet [default: samplesheet.csv]"
    echo -e "  --group_file <path>       Path to group file for sample grouping [header column must have the names 'sample,group' in that order; additional columns may be included, however they will be ignored by the script]"
    echo -e
    "  --group_across_illumina_lanes   Create groups by assuming that files that differ only by a terminal _Lnnn are the same library split across multiple lanes."
    exit 1
fi

if $group_across_illumina_lanes && [[ -n "$group_file" ]]; then
    echo "Provide at most one of --group_file and --group_across_illumina_lanes"
    exit 1
fi

# Display the parameters
echo "Parameters:"
echo "dir_path: $dir_path"
echo "forward_suffix: $forward_suffix"
echo "reverse_suffix: $reverse_suffix"
echo "s3: $s3"
echo "output_path: $output_path"
echo "group_file: $group_file"
echo "group_across_illumina_lanes: $group_across_illumina_lanes"


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

echo "sample,fastq_1,fastq_2" > "$temp_samplesheet"


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
        echo "$sample,${dir_path}${forward_read},${dir_path}${reverse_read}" >> "$temp_samplesheet"
    fi
done

# Check if group file is provided
if [[ -n "$group_file" ]]; then
    # Perform left join with group file
    awk -F',' 'NR==FNR{a[$1]=$2; next} FNR==1{print $0",group"} FNR>1{print $0","(a[$1]?a[$1]:"NA")}' "$group_file" "$temp_samplesheet" > "$output_path"
    echo "CSV file '$output_path' has been created with group information."
elif $group_across_illumina_lanes; then
    cat "$temp_samplesheet" | tr ',' ' ' | \
        while read sample fastq_1 fastq_2; do
            if [[ $sample = "sample" ]]; then
                echo $sample $fastq_1 $fastq_2 "group"
            else
                echo $sample $fastq_1 $fastq_2 \
                     $(echo "$sample" | sed 's/_L[0-9][0-9][0-9]$//')
            fi
        done | tr ' ' ',' > "$output_path"
    echo "CSV file '$output_path' has been created with grouping across illumina lanes."
else
    # If no group file, just use the temporary samplesheet as the final output
    mv "$temp_samplesheet" "$output_path"
    echo "CSV file '$output_path' has been created without group information."
fi

# Remove temporary file if it still exists
[ -f "$temp_samplesheet" ] && rm "$temp_samplesheet"
