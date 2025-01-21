# Use Ubuntu as the base image
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install essential dependencies, R, and system libraries needed for R packages
RUN apt-get update && apt-get install -y \
    wget \
    pigz \
    pbzip2 \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pacman first, then use it to install other packages
RUN R -e "install.packages('pacman', repos='https://cran.rstudio.com/')" && \
    R -e "pacman::p_load(argparse, plotly, ggsci, tidyverse)"

# Set working directory
WORKDIR /opt

# Install Atria
RUN wget https://github.com/cihga39871/Atria/releases/download/v4.1.0/atria-4.1.0-linux-ubuntu22.tar.gz && \
    tar -zxf atria-4.1.0-linux-ubuntu22.tar.gz && \
    rm atria-4.1.0-linux-ubuntu22.tar.gz

# Create symbolic link
RUN chmod +x /opt/atria-4.1.0/bin/atria && \
    ln -s /opt/atria-4.1.0/bin/atria /usr/local/bin/atria

# Set environment variables
ENV PATH="/usr/local/bin:/opt/atria-4.1.0/bin:$PATH"

# Verify installations
RUN which atria && \
    atria --version && \
    Rscript --version && \
    R -e "pacman::p_loaded()"
