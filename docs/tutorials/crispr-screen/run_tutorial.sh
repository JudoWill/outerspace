#!/bin/bash

# CRISPR Screen Analysis Tutorial - Automated Workflow
# This script runs through the complete OUTERSPACE workflow demonstrated in the tutorial
# Uses configuration file approach for consistency and reproducibility

set -e  # Exit on any error

echo "=== OUTERSPACE CRISPR Screen Tutorial ==="
echo "This script demonstrates the configuration file approach for OUTERSPACE workflows."
echo ""

# Create output directories
echo "Creating output directories..."
mkdir -p results/findseq results/collapse results/count results/count_filtered

echo ""
echo "=== Step 1: Extract Sequences (findseq) ==="

# Process control sample (shuffle)
echo "Processing control sample (shuffle)..."
outerspace findseq -c grnaquery.toml \
    -1 data/409-4_S1_L002_R1_001.fastq.gz \
    -2 data/409-4_S1_L002_R2_001.fastq.gz \
    -o results/findseq/shuffle.csv

# Process M1 library sample
echo "Processing M1 library sample..."
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M1_S9_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M1_S9_L001_R2_001.fastq.gz \
    -o results/findseq/M1-lib.csv

# Process M2 library sample
echo "Processing M2 library sample..."
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M2_S12_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M2_S12_L001_R2_001.fastq.gz \
    -o results/findseq/M2-lib.csv

echo ""
echo "=== Step 2: Correct Barcodes (collapse) ==="

# Correct barcodes for all samples using config
echo "Correcting barcodes for all samples using configuration file..."
outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq \
    --output-dir results/collapse

echo ""
echo "=== Step 3: Count Unique Barcodes (count) ==="

# Count barcodes for all samples using config
echo "Counting barcodes for all samples using configuration file..."
outerspace count -c grnaquery.toml \
    --input-dir results/collapse \
    --output-dir results/count

echo ""
echo "=== Step 3b: Count with Allowed List ==="

# Count barcodes using only library protospacers (with config + override)
echo "Counting barcodes with allowed list filter (config + command-line override)..."
outerspace count -c grnaquery.toml \
    --input-dir results/collapse \
    --output-dir results/count_filtered \
    --allowed-list data/library_protospacers.txt

echo ""
echo "=== Step 4: Merge Results ==="

# Merge in wide format using config
echo "Merging results in wide format using configuration file..."
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_wide.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format wide

# Merge in long format using config
echo "Merging results in long format using configuration file..."
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_long.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format long

# Merge filtered results using config
echo "Merging filtered results using configuration file..."
outerspace merge -c grnaquery.toml \
    results/count_filtered/shuffle.csv \
    results/count_filtered/M1-lib.csv \
    results/count_filtered/M2-lib.csv \
    --output-file results/merged_filtered_counts.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format wide

echo ""
echo "=== Step 5: Generate Statistics ==="

# Generate comprehensive statistics using config
echo "Generating statistics for all samples using configuration file..."
outerspace stats -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv

echo ""
echo "=== Demonstrating Command-Line Overrides ==="

# Show how to override config parameters
echo "Demonstrating parameter override: using stricter mismatch tolerance..."
mkdir -p results/collapse_strict

outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq \
    --output-dir results/collapse_strict \
    --mismatches 1 \
    --method cluster

echo "Override complete - results saved to results/collapse_strict/"

echo ""
echo "=== Tutorial Complete! ==="
echo ""
echo "This tutorial demonstrated the configuration file approach with:"
echo "  ✓ Consistent parameter usage across all commands"
echo "  ✓ Easy reproducibility with version-controlled config files"
echo "  ✓ Command-line overrides for testing and customization"
echo ""
echo "Results have been saved to the 'results/' directory:"
echo "  - results/findseq/: Extracted sequences"
echo "  - results/collapse/: Barcode-corrected files (config approach)"
echo "  - results/collapse_strict/: Barcode-corrected files (with override)"
echo "  - results/count/: Barcode counts per protospacer"
echo "  - results/count_filtered/: Filtered barcode counts"
echo "  - results/merged_*.csv: Merged results across samples"
echo ""
echo "Key takeaways:"
echo "  • Configuration files ensure consistency and reproducibility"
echo "  • Command-line parameters can override config settings when needed"
echo "  • Version control your config files for better analysis tracking"
echo ""
echo "You can now explore the results and run additional analyses." 


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.