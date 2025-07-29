FROM condaforge/mambaforge:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI v2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws/

# Install bowtie2 and samtools 
RUN mamba install -c bioconda minimap2=2.28 -y && \
    mamba install -c bioconda samtools=1.21 -y && \
    mamba clean -a

