#!/bin/bash

# Script to safely download databases with file locking
# Usage: download-db.sh <s3_path>
# Example: download-db.sh s3://nao-mgs-index/20250404/output/results/kraken_db

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <s3_path>"
    echo "Example: $0 s3://nao-mgs-index/20250404/output/results/kraken_db"
    exit 1
fi

S3_PATH="$1"

# Extract database name from S3 path (last component)
DB_NAME=$(basename "${S3_PATH}")
LOCAL_PATH="/scratch/${DB_NAME}"

# Create lock file path
LOCK_FILE="/scratch/${DB_NAME}.lock"
COMPLETION_FILE="${LOCAL_PATH}/download_complete.txt"

# Acquire exclusive lock
exec 200>"${LOCK_FILE}"
flock -x 200

# Check if database already downloaded
if [ ! -f "${COMPLETION_FILE}" ]; then
    echo "Downloading ${DB_NAME} from ${S3_PATH} to ${LOCAL_PATH}..."
    
    # Configure AWS S3 settings for optimal transfer
    aws configure set default.s3.max_concurrent_requests 20
    aws configure set default.s3.multipart_threshold 64MB
    aws configure set default.s3.multipart_chunksize 16MB
    
    # Download database
    aws s3 cp --recursive "${S3_PATH}" "${LOCAL_PATH}" && touch "${COMPLETION_FILE}"
    
    echo "Download of ${DB_NAME} completed"
else
    echo "${DB_NAME} already downloaded"
fi

# Release lock
flock -u 200