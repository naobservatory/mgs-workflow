FROM condaforge/mambaforge:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI v2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws/

# Install Kraken2 v2.1.3
# Note that as of 2025-07-29 higher versions segfault with memory-mapping 
RUN mamba install -c bioconda kraken2=2.1.3 -y && \
    mamba clean -a

