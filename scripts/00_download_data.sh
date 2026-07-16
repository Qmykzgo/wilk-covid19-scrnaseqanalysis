#!/bin/bash
# Download script for Wilk et al. 2020 COVID-19 PBMC dataset from CELLxGENE
# Downloads directly to the workspace to avoid stale WSL drive mounts.

DATA_URL="https://datasets.cellxgene.cziscience.com/419da3c2-9141-4654-817f-ee6472df4be3.h5ad"
OUTPUT_FILE="wilk_covid19_pbmc.h5ad"

# Clean up any stale symbolic links
rm -f wilk_covid19_pbmc.h5ad
rm -f results

# Re-create results folder as a local directory
mkdir -p results

echo "Downloading processed Wilk et al. 2020 single-cell dataset (~220MB)..."
if command -v wget &> /dev/null; then
    wget -O "$OUTPUT_FILE" "$DATA_URL"
elif command -v curl &> /dev/null; then
    curl -o "$OUTPUT_FILE" "$DATA_URL"
else
    echo "Error: Neither wget nor curl found. Please install one."
    exit 1
fi

echo "Download completed successfully!"
