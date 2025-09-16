#!/bin/bash

# SIV Barcoding Analysis Tutorial - Automated Workflow
# This script runs through the complete OUTERSPACE workflow for SIV barcode analysis
# Processes BAM files to extract and count viral barcodes across timepoints

set -e  # Exit on any error

echo "=== OUTERSPACE SIV Barcoding Tutorial ==="
echo "This script demonstrates viral barcode extraction and counting from BAM files."
echo ""

# Create output directories
echo "Creating output directories..."
mkdir -p results/findseq results/collapse results/count

echo ""
echo "=== Step 1: Extract Viral Barcodes (findseq) ==="

# Process early timepoint sample
echo "Processing early timepoint sample (MJ21-20230406)..."
outerspace findseq \
    --config sivbarcode.toml \
    -1 data/MJ21-20230406.sifted.bam \
    -o results/findseq/MJ21-20230406.csv

# Process later timepoint sample
echo "Processing later timepoint sample (MJ21-20230711)..."
outerspace findseq \
    --config sivbarcode.toml \
    -1 data/MJ21-20230711.sifted.bam \
    -o results/findseq/MJ21-20230711.csv

echo ""
echo "=== Step 2: Correct Viral Barcodes (collapse) ==="

# Correct barcodes for all samples using config
echo "Correcting viral barcodes for all samples using configuration file..."
outerspace collapse \
    --input-dir results/findseq \
    --output-dir results/collapse \
    --config sivbarcode.toml

echo ""
echo "=== Step 3: Count Unique Viral Barcodes (count) ==="

# Count viral barcodes for all samples
echo "Counting unique viral barcodes for all samples..."
outerspace count \
    --input-dir results/collapse \
    --output-dir results/count \
    --barcode-column viral_barcode_corrected \
    --key-column viral_barcode_corrected

echo ""
echo "=== Step 4: Merge Results ==="

# Merge in wide format for easy comparison
echo "Merging results in wide format for timepoint comparison..."
outerspace merge \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv \
    --output-file results/merged_viral_barcodes.csv \
    --key-column viral_barcode_corrected \
    --count-column viral_barcode_corrected_count \
    --sample-names Early Later \
    --format wide

# Merge in long format for downstream analysis
echo "Merging results in long format for statistical analysis..."
outerspace merge \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv \
    --output-file results/merged_viral_barcodes_long.csv \
    --key-column viral_barcode_corrected \
    --count-column viral_barcode_corrected_count \
    --sample-names Early Later \
    --format long

echo ""
echo "=== Step 5: Generate Statistics ==="

# Generate comprehensive statistics
echo "Generating barcode diversity statistics..."
outerspace stats \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv \
    --key-column viral_barcode_corrected \
    --count-column viral_barcode_corrected_count > results/barcode_statistics.csv

echo ""
echo "=== Analysis Complete! ==="
echo ""
echo "This tutorial demonstrated SIV viral barcode analysis with:"
echo "  ✓ Extraction of viral barcodes from BAM files"
echo "  ✓ Error correction of barcode sequences" 
echo "  ✓ Quantification of viral barcode frequencies"
echo "  ✓ Comparison across timepoints"
echo "  ✓ Statistical analysis of barcode diversity"
echo ""
echo "Results have been saved to the 'results/' directory:"
echo "  - results/findseq/: Extracted viral sequences"
echo "  - results/collapse/: Error-corrected viral barcodes"
echo "  - results/count/: Viral barcode counts per sample"
echo "  - results/merged_viral_barcodes.csv: Combined results (wide format)"
echo "  - results/merged_viral_barcodes_long.csv: Combined results (long format)"
echo "  - results/barcode_statistics.csv: Diversity statistics"
echo "  - results/barcode_summary.csv: Summary table of barcode counts"
echo ""

echo ""
echo "Next steps for analysis:"
echo "  1. Examine barcode_statistics.csv for diversity metrics"
echo "  2. Analyze temporal dynamics in merged_viral_barcodes.csv"
echo "  3. Identify dominant viral variants at each timepoint"
echo "  4. Correlate barcode patterns with biological outcomes"
echo ""
echo "Tutorial complete! Explore the results directory for detailed outputs."


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
