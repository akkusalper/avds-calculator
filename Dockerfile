FROM python:3.11-slim

LABEL maintainer="A. Akkus, O. Ozdemir - Acibadem Mehmet Ali Aydinlar University"
LABEL description="Allele-Specific Variant Detection Score Calculator"
LABEL version="1.0"
LABEL organization="Department of Translational Medicine, Acibadem Mehmet Ali Aydinlar University, Istanbul, Turkey"

# Install system dependencies for bioinformatics tools
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first (for Docker layer caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY avds_calculator.py .
COPY avds_pipeline.py .
COPY README.md .

# Create data directories
RUN mkdir -p /data/input /data/output

# Set Python to run in unbuffered mode (better for logging)
ENV PYTHONUNBUFFERED=1

# Set entrypoint to the pipeline script
ENTRYPOINT ["python", "avds_pipeline.py"]

# Default command shows help
CMD ["--help"]
