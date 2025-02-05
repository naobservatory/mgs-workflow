FROM debian:bullseye-slim

# Install Python and dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    gcc \
    libffi-dev \
    libssl-dev \
    zlib1g \
    libbz2-dev \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir pandas biopython

# Set the working directory
WORKDIR /data
