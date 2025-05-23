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

## Approach 1: Using the Pipeline Command

The pipeline command automates all steps in one command. This is the simplest approach:

```bash
# Create output directories
mkdir -p tmp/

# Run the pipeline
outerspace pipeline tests/configs/grnaquery.cfg \
    --input-dir tests/data \
    --output-dir tmp/results \
    --barcode-columns UMI_5prime,UMI_3prime \
    --key-column protospacer \
    --mismatches 2 \
    --method directional \
    --metrics
```

This will:
1. Process all FASTQ files in `tests/data`
2. Extract sequences using patterns in `grnaquery.cfg`
3. Correct barcodes using UMI-tools clustering
4. Count unique barcodes per protospacer
5. Generate metrics files for quality control

The results will be organized in the `tmp/results` directory:
- `extracted/`: Raw sequence extractions
- `collapsed/`: Barcode-corrected files
- `counted/`: Final barcode counts
- `collapse_metrics.yaml`: Metrics from barcode correction
- `count_metrics.yaml`: Metrics from counting

## Approach 2: Using Individual Commands

For more control over each step, you can run the commands individually:

### Step 1: Extract Sequences

First, extract sequences from each sample's FASTQ files:

```bash
# Create output directory
mkdir -p tmp/results/extracted

# Process control sample
outerspace findseq tests/configs/grnaquery.cfg \
    -1 tests/data/409-4_S1_L002_R1_001.fastq.gz \
    -2 tests/data/409-4_S1_L002_R2_001.fastq.gz \
    -o tmp/results/extracted/409-4_S1_L002_R1_001p.csv

# Process M1 sample
outerspace findseq tests/configs/grnaquery.cfg \
    -1 tests/data/2-G1L9-M1_S9_L001_R1_001.fastq.gz \
    -2 tests/data/2-G1L9-M1_S9_L001_R2_001.fastq.gz \
    -o tmp/results/extracted/2-G1L9-M1_S9_L001_R1_001p.csv

# Process M2 sample
outerspace findseq tests/configs/grnaquery.cfg \
    -1 tests/data/2-G1L9-M2_S12_L001_R1_001.fastq.gz \
    -2 tests/data/2-G1L9-M2_S12_L001_R2_001.fastq.gz \
    -o tmp/results/extracted/2-G1L9-M2_S12_L001_R1_001p.csv
```

### Step 2: Correct Barcodes

Next, correct sequencing errors in the barcodes:

```bash
# Create output directory
mkdir -p tmp/results/collapsed

# Run collapse command on all extracted files
outerspace collapse \
    --input-dir tmp/results/extracted \
    --output-dir tmp/results/collapsed \
    --columns UMI_5prime,UMI_3prime \
    --mismatches 2 \
    --method directional \
    --metrics tmp/results/collapsed/collapse_metrics.yaml
```

### Step 3: Count Barcodes

Finally, count unique barcodes per protospacer:

```bash
# Create output directory
mkdir -p tmp/results/counted

# Run count command on all collapsed files
outerspace count \
    --input-dir tmp/results/collapsed \
    --output-dir tmp/results/counted \
    --barcode-column UMI_5prime_UMI_3prime_corrected \
    --key-column protospacer \
    --metrics tmp/results/counted/count_metrics.yaml
```

### Optional: Calculate Gini Coefficients

You can calculate Gini coefficients to assess the distribution of barcodes:

```bash
# Calculate Gini for control sample
outerspace gini tmp/results/counted/409-4_S1_L002_R1_001p.csv \
    --column counts

# Calculate Gini for M1 sample
outerspace gini tmp/results/counted/2-G1L9-M1_S9_L001_R1_001p.csv \
    --column counts

# Calculate Gini for M2 sample
outerspace gini tmp/results/counted/2-G1L9-M2_S12_L001_R1_001p.csv \
    --column counts
```

### Optional: Visualize Results

Create visualizations of the barcode distributions:

```bash
# Create output directory
mkdir -p tmp/results/plots

# Generate plots
outerspace visualize \
    tmp/results/counted \
    tmp/results/plots \
    --title-prefix "CRISPR Screen" \
    --log-scale
```

## Comparing the Approaches

The pipeline command is simpler and ensures consistent processing across all files. However, using individual commands gives you more control and allows you to:

1. Inspect intermediate results between steps
2. Modify parameters for specific samples
3. Rerun individual steps if needed
4. Process files in parallel
5. Add custom processing steps

Choose the approach that best fits your needs. For routine processing, the pipeline command is recommended. For more complex analyses or troubleshooting, individual commands provide more flexibility.

## Next Steps

After processing, you can:
1. Compare barcode distributions between samples
2. Identify enriched or depleted protospacers
3. Analyze the metrics files for quality control
4. Generate publication-quality visualizations
5. Perform statistical analysis on the results 