# Complete Spatial Transcriptomics Pipeline
FROM python:3.11-slim

WORKDIR /workspace

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies including Java for Nextflow
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    curl \
    wget \
    git \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/default-java

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

# Copy requirements and install Python packages
COPY docker/requirements.txt /workspace/
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Create directories
RUN mkdir -p /workspace/{data,results,logs,modules,workflows,notebooks}

# Copy pipeline code (commented out for now - add when you have the files)
# COPY modules/ /workspace/modules/
# COPY workflows/ /workspace/workflows/
# COPY notebooks/ /workspace/notebooks/

# Set Python path
ENV PYTHONPATH="/workspace"

# Default command
CMD ["/bin/bash"]