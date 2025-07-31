#!/bin/bash

# Script to safely download databases with file locking
# Usage: download-db.sh <s3_path> [timeout_seconds]
# Example: download-db.sh s3://nao-mgs-index/20250404/output/results/kraken_db 1200

set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <s3_path> [timeout_seconds]"
    echo "Example: $0 s3://nao-mgs-index/20250404/output/results/kraken_db 1200"
    exit 1
fi

# Set timeout (use provided value or no timeout if not specified)
TIMEOUT_SECONDS=${2:-""}

# Normalize S3 path (replace consecutive slashes with a single slash, then add double slash in s3:// protocol back in)
S3_PATH=$(echo "$1" | sed -e 's|///*|/|g' -e 's|^s3:/|s3://|')

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
if [ -n "$TIMEOUT_SECONDS" ]; then
    flock -x -w "$TIMEOUT_SECONDS" 200 || { echo "Timed out waiting for lock after $TIMEOUT_SECONDS seconds"; exit 1; }
else
    flock -x 200
fi

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
