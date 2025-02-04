# Use Amazon Linux as the base image
FROM amazonlinux:latest

# Install Docker and other necessary packages
RUN yum update -y && \
    yum install -y docker unzip zip tar bash findutils wget git && \
    yum clean all

# Install SDKMAN and Java 17
RUN curl -s "https://get.sdkman.io" | bash && \
    bash -c "source /root/.sdkman/bin/sdkman-init.sh && sdk install java 17.0.10-tem"

# Persist Java and Nextflow in the PATH
ENV SDKMAN_DIR="/root/.sdkman"
ENV PATH="${SDKMAN_DIR}/candidates/java/current/bin:/root/.local/bin:$PATH"

# Verify Java installation
RUN java -version

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mkdir -p /root/.local/bin/ && \
    mv nextflow /root/.local/bin && \
    nextflow info

# Install nf-test
RUN wget -qO- https://get.nf-test.com | bash && \
    mv nf-test /usr/local/bin/

# Ensure the user starts in the home directory
WORKDIR /root

# Download our repo
RUN git clone https://github.com/naobservatory/mgs-workflow

# Set the default command
CMD ["/bin/bash"]
