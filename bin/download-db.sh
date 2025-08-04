#!/bin/bash

# Script to safely download/copy databases with file locking
# Usage: download-db.sh <source_path> [timeout_seconds]
# Examples: 
#   download-db.sh s3://nao-mgs-index/20250404/output/results/kraken_db 1200
#   download-db.sh /path/to/local/kraken_db 1200

set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <source_path> [timeout_seconds]"
    echo "Examples:"
    echo "  $0 s3://nao-mgs-index/20250404/output/results/kraken_db 1200"
    echo "  $0 /path/to/local/kraken_db 1200"
    exit 1
fi

# Set timeout (use provided value or no timeout if not specified)
TIMEOUT_SECONDS=${2:-""}

# Handle both S3 and local paths
SOURCE_PATH="$1"

# Normalize path (replace consecutive slashes with a single slash)
SOURCE_PATH=$(echo "$SOURCE_PATH" | sed -e 's|///*|/|g')

# Check if this is an S3 path and restore s3:// protocol if needed (above normalization will have make it start with s3:/)
if [[ "$SOURCE_PATH" =~ ^s3: ]]; then
    IS_S3=true
    SOURCE_PATH="s3://${SOURCE_PATH#s3:/}"
else
    IS_S3=false
fi

# Extract database name from source path (last component)
DB_NAME=$(basename "${SOURCE_PATH}")
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

# Check if database already downloaded/copied
if [ ! -f "${COMPLETION_FILE}" ]; then
    if [ "$IS_S3" = true ]; then
        echo "Downloading ${DB_NAME} from ${SOURCE_PATH} to ${LOCAL_PATH}..."
        
        # Configure AWS S3 settings for optimal transfer
        aws configure set default.s3.max_concurrent_requests 20
        aws configure set default.s3.multipart_threshold 64MB
        aws configure set default.s3.multipart_chunksize 16MB
        
        # Download database
        aws s3 sync "${SOURCE_PATH}" "${LOCAL_PATH}" && touch "${COMPLETION_FILE}"
    else
        echo "Copying ${DB_NAME} from ${SOURCE_PATH} to ${LOCAL_PATH}..."
        
        # Check if source exists
        if [ ! -e "${SOURCE_PATH}" ]; then
            echo "Error: Source path ${SOURCE_PATH} does not exist"
            exit 1
        fi
        
        # Copy database recursively
        cp -r "${SOURCE_PATH}/." "${LOCAL_PATH}" && touch "${COMPLETION_FILE}"
    fi
    
    echo "Transfer of ${DB_NAME} completed"
else
    echo "${DB_NAME} already available"
fi
