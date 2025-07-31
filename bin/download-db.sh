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

# Normalize S3 path (remove double slashes, preserve s3:// protocol)
temp="${1#s3://}" # Remove s3:// prefix
temp="${temp//\/\/*/\/}"  # Replace multiple slashes with single
S3_PATH="s3://${temp#/}"  # Remove leading slash and add protocol

# Extract database name from S3 path (last component)
DB_NAME=$(basename "${S3_PATH}")
LOCAL_PATH="/scratch/${DB_NAME}"

mkdir -p /scratch

# Create lock file path
LOCK_FILE="/scratch/${DB_NAME}.lock"
COMPLETION_FILE="${LOCAL_PATH}/download_complete.txt"

# We want to make sure lock is released even if the script exits unexpectedly
# Note that the trap will also handle release on successful completion
cleanup() {
  flock -u 200 || true
  exec 200>&- 
}
trap cleanup EXIT ERR INT TERM


# Acquire exclusive lock
exec 200>"${LOCK_FILE}"
flock -x -w 1200 200 || { echo "Timed out waiting for lock"; exit 1; }

# Check if database already downloaded
if [ ! -f "${COMPLETION_FILE}" ]; then
    echo "Downloading ${DB_NAME} from ${S3_PATH} to ${LOCAL_PATH}..."
    
    # Configure AWS S3 settings for optimal transfer
    aws configure set default.s3.max_concurrent_requests 20
    aws configure set default.s3.multipart_threshold 64MB
    aws configure set default.s3.multipart_chunksize 16MB
    
    # Download database
    aws s3 sync "${S3_PATH}" "${LOCAL_PATH}" && touch "${COMPLETION_FILE}"
    
    echo "Download of ${DB_NAME} completed"
else
    echo "${DB_NAME} already downloaded"
fi
