# SIV Barcoding Analysis Tutorial

This tutorial demonstrates how to analyze SIV barcoding data using OUTERSPACE's individual commands. We'll process BAM files containing viral barcode sequences to extract and count unique viral barcodes across different samples.

## Overview

We'll analyze two samples from a longitudinal SIV study:
- `MJ21-20230406.sifted.bam`: Early timepoint sample
- `MJ21-20230711.sifted.bam`: Later timepoint sample

Each BAM file contains reads that may include:
- **viral_barcode**: 34-nucleotide unique barcode within the viral genome
- **left_umi**: 18-nucleotide UMI sequence on the left flank
- **right_umi**: 18-nucleotide UMI sequence on the right flank
- **LuttegePolP**: Specific viral sequence marker

## Prerequisites

1. Clone the OUTERSPACE repository:
```bash
git clone https://github.com/DamLabResources/outerspace.git
cd outerspace
```

2. Install OUTERSPACE:
```bash
pip install -e .
```

3. Navigate to this tutorial directory:
```bash
cd docs/tutorials/siv-barcoding
```

## Data Files

This tutorial includes the following data files in the `data/` directory:
- `MJ21-20230406.sifted.bam`: Early timepoint BAM file
- `MJ21-20230711.sifted.bam`: Later timepoint BAM file

## Understanding the Configuration File

OUTERSPACE uses TOML configuration files to define patterns and command parameters. Let's break down the `sivbarcode.toml` file:

### Global Patterns Section

```toml
[[patterns]]
name = "viral_barcode"
reg_expr = "(cccctccaggactagcataa){i<=1,d<=1,s<=2}(?P<viral_barcode>.{34})(atggaagaaagacctccaga){i<=1,d<=1,s<=2}"
read = "R1"
orientation = "both"
multiple = "first"
```

This defines a pattern to extract the viral barcode:
- **`name`**: Identifier for the pattern, becomes a column name in output
- **`reg_expr`**: Regular expression with named capture group `(?P<viral_barcode>.{34})`
  - `.{34}` captures exactly 34 nucleotides for the viral barcode
  - Flanked by constant sequences with fuzzy matching allowing insertions, deletions, and substitutions
  - `{i<=1,d<=1,s<=2}` allows up to 1 insertion, 1 deletion, and 2 substitutions
- **`read`**: Which read to search (R1 for single-end BAM files)
- **`orientation`**: Search direction (both forward and reverse)
- **`multiple`**: How to handle multiple matches (first match only)

```toml
[[patterns]]
name = "left_umi"
reg_expr = "(CAAGCAGAAGACGGCATACGAGAT){i<=1,d<=1,s<=2}(?<left_umi>.{18})(atatacttagaaaaggaagaaggcatcataccaga){i<=1,d<=1,s<=2}"
read = "R1"
orientation = "both"
multiple = "first"
```

The left UMI pattern extracts an 18-nucleotide barcode:
- Captures exactly 18 nucleotides between constant flanking sequences
- Uses fuzzy matching for the flanking regions

```toml
[[patterns]]
name = "right_umi"
reg_expr = "(gaagaactccgaaaaaggctaaggc){i<=1,d<=1,s<=2}(?<right_umi>.{18})(GATCTCGGTGGTCGCCGTATCATT){i<=1,d<=1,s<=2}"
read = "R1"
orientation = "both"
multiple = "first"
```

The right UMI pattern extracts another 18-nucleotide barcode from a different region.

```toml
[[patterns]]
name = "LuttegePolP"
reg_expr = "(?<LuttegePolP>tggggtaccatacaatccac){i<=1,d<=1,s<=2}"
read = "R1"
orientation = "both"
multiple = "first"
```

The LuttegePolP pattern identifies the 'probe sequence' for a common marker of viral intactness.

### Command-Specific Configurations

```toml
[findseq]
use_all_patterns = true
matches_only = true
threads = 8
skip_unmapped = true
```

- **`use_all_patterns`**: Extract all defined patterns from each read
- **`matches_only`**: Only output reads where at least one pattern is found
- **`threads`**: Use 8 threads for parallel processing
- **`skip_unmapped`**: Skip unmapped reads in BAM files

```toml
[collapse]
columns = 'viral_barcode'
mismatches = 2
```

- **`columns`**: Only apply barcode correction to the viral_barcode column
- **`mismatches`**: Maximum edit distance for clustering similar barcodes

## Tutorial Workflow (Configuration File Approach)

This is the **recommended approach** as it ensures consistency and reproducibility across all steps.

### Step 1: Extract Sequences (findseq)

Extract viral barcodes and other sequences from BAM files using the configured patterns:

```bash
# Create output directories
mkdir -p results/findseq results/collapse results/count

# Process early timepoint sample
outerspace findseq -c sivbarcode.toml \
    -1 data/MJ21-20230406.sifted.bam \
    -o results/findseq/MJ21-20230406.csv

# Process later timepoint sample  
outerspace findseq -c sivbarcode.toml \
    -1 data/MJ21-20230711.sifted.bam \
    -o results/findseq/MJ21-20230711.csv
```

**Expected output**: CSV files containing extracted viral barcodes and UMI sequences for each sample.

### Step 2: Correct Barcodes (collapse)

Correct sequencing errors in the viral barcodes using the configuration file settings:

```bash
# Correct barcodes for all samples using config
outerspace collapse -c sivbarcode.toml \
    --input-dir results/findseq \
    --output-dir results/collapse
```

**What this does**:
- Uses `columns = 'viral_barcode'` from config to only correct viral barcodes
- Applies `mismatches = 2` from config for clustering
- Creates corrected barcode column: `viral_barcode_corrected`

### Step 3: Count Unique Barcodes (count)

Count unique viral barcodes per sample:

```bash
# Count viral barcodes for all samples (uses config settings)
outerspace count -c sivbarcode.toml \
    --input-dir results/collapse \
    --output-dir results/count
```

**Expected output**: CSV files showing the count of each unique viral barcode in each sample.

### Step 4: Merge Results

Combine all sample results for comparative analysis:

```bash
# Merge in wide format for easy comparison (uses config settings)
outerspace merge -c sivbarcode.toml \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv \
    --output-file results/merged_viral_barcodes.csv \
    --sample-names Early Later \
    --format wide

# Merge in long format for downstream analysis (uses config settings)
outerspace merge -c sivbarcode.toml \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv \
    --output-file results/merged_viral_barcodes_long.csv \
    --sample-names Early Later \
    --format long
```

### Step 5: Generate Statistics

Calculate summary statistics for barcode diversity:

```bash
# Generate comprehensive statistics (uses config settings)
outerspace stats -c sivbarcode.toml \
    results/count/MJ21-20230406.csv \
    results/count/MJ21-20230711.csv > results/barcode_statistics.csv
```

## Single File Processing

You can also process files individually for more control:

```bash
# Process a single sample through the entire workflow
outerspace findseq -c sivbarcode.toml \
    -1 data/MJ21-20230406.sifted.bam \
    -o results/single_early.csv

outerspace collapse -c sivbarcode.toml \
    --input-file results/single_early.csv \
    --output-file results/single_early_collapsed.csv

outerspace count -c sivbarcode.toml \
    --input-file results/single_early_collapsed.csv \
    --output-file results/single_early_counts.csv
```

## Understanding the Results

### findseq Output
- Each row represents one read from the BAM file
- Columns contain extracted sequences (viral_barcode, left_umi, right_umi, LuttegePolP)
- Only reads with at least one pattern match are included (due to `matches_only = true`)

### collapse Output
- Adds corrected barcode column (`viral_barcode_corrected`)
- Groups similar viral barcodes to reduce sequencing errors
- Preserves original barcode information for comparison

### count Output
- Each row represents one unique viral barcode
- Shows count of occurrences per barcode
- Higher counts indicate more abundant viral variants

### merge Output
- **Wide format**: One row per viral barcode, one column per timepoint
- **Long format**: One row per barcode-timepoint combination
- Facilitates comparison of barcode frequencies between timepoints

## Quality Control and Interpretation

Monitor these metrics throughout the analysis:

1. **Extraction efficiency**: Percentage of reads with viral patterns found
2. **Barcode diversity**: Number of unique viral barcodes per sample
3. **Temporal dynamics**: Changes in barcode frequencies between timepoints
4. **Barcode correction**: Reduction in unique barcodes after clustering

### Key Questions for SIV Barcode Analysis:

1. **Viral diversity**: How many unique viral variants are present?
2. **Persistence**: Which viral barcodes persist across timepoints?
3. **Expansion/contraction**: Which variants increase or decrease over time?
4. **Bottlenecks**: Are there evidence of population bottlenecks?

## Command-Line Overrides

You can override specific config settings via command-line parameters when needed for testing or customization.

## Best Practices

1. **Use configuration files** for reproducible analyses
2. **Version control your config files** to track parameter changes
3. **Validate patterns** with a subset of data before full processing
4. **Compare timepoints** to identify viral population dynamics
5. **Document analysis decisions** for collaborative research

## Next Steps

After completing this tutorial, you can:
1. Analyze barcode diversity trends over time
2. Identify dominant viral variants at each timepoint
3. Calculate diversity indices (Shannon, Simpson) for viral populations
4. Correlate barcode dynamics with clinical outcomes
5. Perform statistical analysis to identify significant changes

## Troubleshooting

**Low extraction rates**: Check BAM file quality and pattern configurations
**High barcode diversity**: Consider adjusting mismatch tolerance in collapse step
**Memory issues**: Process samples individually or use smaller chunks
**Missing patterns**: Verify viral sequences match expected barcode contexts
**Parameter conflicts**: Command-line arguments always override config file settings

Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
