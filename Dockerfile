# Use Python 3.11 slim image for smaller size
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy requirements file
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host files.pythonhosted.org -r requirements.txt

# Copy the application code
COPY demux_barcodes.py .

# Create output directory
RUN mkdir -p /data/output

# Set the working directory for data
WORKDIR /data

# Set entrypoint to python script
ENTRYPOINT ["python", "/app/demux_barcodes.py"]

# Default command (shows help if no arguments provided)
CMD ["--help"]
