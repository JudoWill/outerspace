# OUTERSPACE Walkthrough

This walkthrough demonstrates how to use OUTERSPACE to process CRISPR screen data. We'll show two approaches:
1. Using the `pipeline` command to process everything in one step
2. Using individual commands for more control over each step

## Data Overview

We'll be processing three samples from a CRISPR screen:
- `409-4_S1_L002`: A single gRNA
- `2-G1L9-M1_S9_L001`: A first generation library (M1)
- `2-G1L9-M2_S12_L001`: A second generation libary (M2)

Each sample has paired-end FASTQ files (R1 and R2) containing the sequencing reads.

## Configuration File

Create a configuration file `config.toml` with the new global patterns structure:

```toml
# Global patterns
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}"
read = "R1"
orientation = "forward"
multiple = "first"

[[patterns]]
name = "protospacer"
reg_expr = "(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})"
read = "R1"
orientation = "forward"
multiple = "first"

[[patterns]]
name = "UMI_3prime"
reg_expr = "(?P<UMI_3prime>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}"
read = "R2"
orientation = "forward"
multiple = "first"

[findseq]
pattern_names = ["UMI_5prime", "protospacer", "UMI_3prime"]
```

## Approach 1: Using the Pipeline Command

The pipeline command automates all steps in one command. This is the simplest approach:

```bash
# Create output directories
mkdir -p tmp/

# Run the pipeline
outerspace pipeline config.toml snakemake_config.yaml \
    --snakemake-args="--cores 4"
```

This will:
1. Process all FASTQ files specified in the Snakemake configuration
2. Extract sequences using patterns in `config.toml`
3. Correct barcodes using UMI-tools clustering
4. Count unique barcodes per protospacer
5. Generate metrics files for quality control

The results will be organized in the output directory structure defined in your Snakemake configuration:
- `findseq/`: Raw sequence extractions
- `collapse/`: Barcode-corrected files
- `count/`: Final barcode counts per sample
- `merge/`: Merged results across samples
- `stats/`: Statistical summaries and metrics

See the [Snakemake](docs/snakemake.md) documentation for more details about how to configure the pipeline.

## Approach 2: Using Individual Commands

For more control over each step, you can run the commands individually:

### Step 1: Extract Sequences

First, extract sequences from each sample's FASTQ files:

```bash
# Create output directory
mkdir -p tmp/results/findseq

# Process control sample
outerspace findseq -c config.toml \
    -1 tests/data/409-4_S1_L002_R1_001.fastq.gz \
    -2 tests/data/409-4_S1_L002_R2_001.fastq.gz \
    -o tmp/results/findseq/409-4_S1_L002.csv

# Process M1 sample
outerspace findseq -c config.toml \
    -1 tests/data/2-G1L9-M1_S9_L001_R1_001.fastq.gz \
    -2 tests/data/2-G1L9-M1_S9_L001_R2_001.fastq.gz \
    -o tmp/results/findseq/2-G1L9-M1_S9_L001.csv

# Process M2 sample
outerspace findseq -c config.toml \
    -1 tests/data/2-G1L9-M2_S12_L001_R1_001.fastq.gz \
    -2 tests/data/2-G1L9-M2_S12_L001_R2_001.fastq.gz \
    -o tmp/results/findseq/2-G1L9-M2_S12_L001.csv
```

### Step 2: Correct Barcodes

Next, correct sequencing errors in the barcodes:

```bash
# Create output directory
mkdir -p tmp/results/collapse

# Run collapse command on all extracted files (uses config settings)
outerspace collapse -c config.toml \
    --input-dir tmp/results/findseq \
    --output-dir tmp/results/collapse \
    --metrics tmp/results/collapse/collapse_metrics.yaml
```

### Step 3: Count Barcodes

Finally, count unique barcodes per protospacer:

```bash
# Create output directory
mkdir -p tmp/results/count

# Run count command on all collapsed files (uses config settings)
outerspace count -c config.toml \
    --input-dir tmp/results/collapse \
    --output-dir tmp/results/count
```

### Step 4: Merge Results (Optional)

You can merge all count files into a single file for downstream analysis:

```bash
# Create output directory
mkdir -p tmp/results/merge

# Merge all count files (uses config settings)
outerspace merge -c config.toml \
    tmp/results/count/*.csv \
    --output-file tmp/results/merge/merged_counts.csv \
    --sample-names 409-4 M1 M2 \
    --format wide
```

### Optional: Calculate Statistics

You can calculate comprehensive statistics including Gini coefficients, Shannon diversity, and other metrics:

```bash
# Calculate statistics for all samples (uses config settings)
outerspace stats -c config.toml \
    tmp/results/count/*.csv > tmp/results/count/statistics.csv
```

### Optional: Visualize Results

Create visualizations of the barcode distributions:

```bash
# Create output directory
mkdir -p tmp/results/plots

# Generate plots (uses config settings for styling)
outerspace visualize -c config.toml \
    tmp/results/count \
    tmp/results/plots \
    --title-prefix "CRISPR Screen" \
    --log-scale
```

## Next Steps

After processing, you can:
1. Compare barcode distributions between samples
2. Identify enriched or depleted protospacers
3. Analyze the metrics files for quality control
4. Generate publication-quality visualizations
5. Perform statistical analysis on the results 


Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
