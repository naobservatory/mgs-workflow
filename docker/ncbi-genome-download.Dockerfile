FROM debian:bullseye-slim

# Install Python and dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    procps \
    gcc \
    libffi-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install ncbi-genome-download
RUN pip3 install --no-cache-dir ncbi-genome-download

# Ensure ncbi-genome-download is accessible as a shell command
ENV PATH="/usr/local/bin:$PATH"