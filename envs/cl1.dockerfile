# Use the latest Debian base image
FROM debian:bookworm-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install system dependencies
RUN apt-get update && apt-get install -y \
    default-jre \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Download and install ClusterOne 1.2
#RUN wget https://github.com/ntamas/cl1/archive/refs/tags/1.2.tar.gz
#RUN tar -xvzf 1.2.tar.gz
RUN wget https://github.com/ntamas/cl1/releases/download/1.2/cluster_one-1.2.jar