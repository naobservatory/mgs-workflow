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

# Download any plugins as running in parallel doesn't guarantee that they install properly
nf-test update-plugins
sleep 10

# Create sequence of shard numbers (1 to total_threads)
for shard in $(seq 1 "$total_threads"); do
    (
        echo "Running shard ${shard}/${total_threads}..."
        nf-test test --shard "${shard}/${total_threads}" "$@"
    ) &
    pids+=($!) # Store process ID
done

# Wait for all background processes
for pid in "${pids[@]}"; do
    wait "$pid"
    statuses+=($?)
done

# Check exit statuses
all_passed=true
for status in "${statuses[@]}"; do
    if [ "$status" -ne 0 ]; then
        all_passed=false
        echo "One or more test shards failed with exit status: $status"
    fi
done

if $all_passed; then
    echo "All test shards completed successfully"
    exit 0
else
    exit 1
fi
