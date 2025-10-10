# Use the latest Debian base image
FROM debian:bookworm-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PIP_BREAK_SYSTEM_PACKAGES=1

# Update package list and install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    libhdf5-dev \
    libfreetype6-dev \
    libpng-dev \
    pkg-config \
    git \
    wget \
    curl \
    zip \
    && rm -rf /var/lib/apt/lists/*

# Create a symbolic link for python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install python packages
RUN python3 -m pip install \
    h5py==3.12.1 \
    openpyxl==3.1.5 \
    numpy==2.1.3 \
    pandas==2.2.3 \
    dask==2025.2.0 \
    pyarrow==19.0.0 \
    scipy==1.15.1 \
    matplotlib==3.10.0 \
    matplotlib-venn==1.1.1 \
    seaborn==0.13.2 \
    scikit-learn==1.6.1 \
    shap==0.46.0 \
    ccc-coef==0.2.2 \
    xarray==2025.7.1 \
    missingno==0.5.2 \
    xgboost==2.1.4 \
    git+https://github.com/gregversteeg/NPEET.git
