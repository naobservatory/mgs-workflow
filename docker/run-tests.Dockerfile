# Use Amazon Linux as the base image
FROM amazonlinux:latest

# Install Docker and other necessary packages
RUN yum update -y && \
    yum install -y docker && \
    yum install -y unzip

RUN yum install -y zip
RUN yum install -y tar
RUN yum install -y bash
RUN yum install -y findutils

# Install SDKMAN and Java 17
RUN curl -s "https://get.sdkman.io" | bash && \
  bash -c "source /root/.sdkman/bin/sdkman-init.sh && sdk install java 17.0.10-tem"

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mkdir -p /root/.local/bin/ && \
    mv nextflow /root/.local/bin && \
    /root/.local/bin/nextflow info

# Install nf-test
RUN wget -qO- https://get.nf-test.com | bash && \
    mv nf-test /usr/local/bin/

# Add SDKMAN to the PATH for interactive use
ENV PATH="/root/.sdkman/candidates/java/current/bin:$PATH"

# Set the default command
CMD ["/bin/bash"]
