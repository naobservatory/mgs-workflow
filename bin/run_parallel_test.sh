#!/bin/bash

# Usage: ./parallel_nf_test.sh <total_threads> [additional_args]
# Example: ./parallel_nf_test.sh 4 --verbose --some-flag

# Function to print usage
print_usage() {
    echo "Usage: $0 <total_threads> [additional_args]"
    echo "Example: $0 4 --verbose --some-flag"
    exit 1
}

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    print_usage
fi

total_threads=$1
shift  # Shift arguments so $@ contains additional args

# Validate input is a positive integer
if ! [[ "$total_threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Please provide a positive integer for total threads"
    print_usage
fi

# Store exit statuses
declare -a pids
declare -a statuses

temp_dir=$(mktemp -d)
error_log="./test_logs.txt"
touch "$error_log"

# Download any plugins as running in parallel doesn't guarantee that they install properly
nf-test update-plugins
sleep 2

# Create sequence of shard numbers (1 to total_threads)
for shard in $(seq 1 "$total_threads"); do
    (
        echo "Running shard ${shard}/${total_threads}..."
        nf-test test --shard "${shard}/${total_threads}" "$@" &> "$temp_dir/error_${shard}.log"
        exit_code=$?
        if [ $exit_code -ne 0 ]; then
            echo "Shard ${shard} failed with exit code $exit_code" >> "$error_log"
            cat "$temp_dir/error_${shard}.log" >> "$error_log"
            echo "---" >> "$error_log"
        fi
        exit $exit_code
    ) &
    pids+=($!) # Store process ID
done

# Wait for all background processes
for pid in "${pids[@]}"; do
    wait "$pid"
    statuses+=($?)
    
    # Collect failed shards
    if [ ${statuses[-1]} -ne 0 ]; then
        failed_shards=true
    fi
done
# Remove temp dir
rm -rf "$temp_dir"

# Check exit statuses
if [ "$failed_shards" == "true" ]; then
    echo "Some test shards failed. Check the test log at: $error_log"
    grep "FAILED" "$error_log"
    exit 1
else
    echo "All test shards completed successfully"
    exit 0
fi
