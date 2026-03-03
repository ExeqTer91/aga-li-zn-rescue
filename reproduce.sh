#!/bin/bash

set -e

# Data directory setup
mkdir -p data/raw data/processed results/tables results/figures results/logs

echo "Step 0: Download GEO"
# Note: Provide your GSE ID here
# python src/00_download_geo.py --gse GSEXXXXX

echo "Step 1: Preprocess & QC"
# python src/01_preprocess_qc.py

echo "Step 2: Differential Expression Analysis"
# python src/02_deg_analysis.py

echo "Pipeline (Part 1) finished successfully!"
