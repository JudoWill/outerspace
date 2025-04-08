# Usage

## CRISPR Screens

OUTERSPACE (Optimized Utilities for Tracking Enrichment in Screens through Precise Analysis of CRISPR Experiments) is designed to analyze data from CRISPR screens. These screens typically involve:

1. Creating a library of guide RNAs (gRNAs) targeting genes of interest
2. Delivery of the gRNA library into cells by lentiviral infection
3. Selection or screening process
4. Sequencing of the gRNA sequences before and after selection
5. Analysis of gRNA abundance changes

Check out this [AddGene writeup](https://www.addgene.org/guides/pooled-libraries/) for more details.

The proposed OUTERSPACE workflow helps process and analyze this data through several steps:

1. **Sequence Extraction** (`find`): Extracts gRNA sequences and associated barcodes from FASTQ files using configurable search patterns.

2. **Barcode Correction** (`collapse`): Corrects sequencing errors in barcodes using UMI-tools algorithms to improve accuracy.

3. **Counting** (`count`): Quantifies the frequency of each gRNA-barcode combination.

4. **Distribution Analysis** (`gini`): Calculates Gini coefficients to measure the inequality of barcode distributions. This is useful for diagnostic purposes.

5. **Difference** (`diff`): Calculates the difference between a before/after comparison. Identifying successful gRNAs.

6. **Visualization** (`visualize`): Creates plots to help interpret the results.

his integrated toolset helps researchers accurately measure changes in gRNA abundance between conditions, identify hits from their screens, and ensure data quality through multiple analysis steps.

### Usage Story

```bash
# Extract sequences
outerspace find config.toml -1 ctrl_r1.fastq -2 ctrl_r2.fastq -o ctrl_output.csv
outerspace find config.toml -1 exp_r1.fastq -2 exp_r2.fastq -o exp_output.csv

# Correct barcodes
outerspace collapse --columns umi3,umi5 --mismatches 2 ctrl_output.csv ctrl_output_collapsed.csv
outerspace collapse --columns umi3,umi5 --mismatches 2 exp_output.csv exp_output_collapsed.csv

# Count barcodes
outerspace count --barcode-column umi3_umi5_corrected --key-column protospacer ctrl_output_collapsed.csv ctrl_counts.csv
outerspace count --barcode-column umi3_umi5_corrected --key-column protospacer exp_output_collapsed.csv exp_counts.csv 

# Calculate Gini coefficient
outerspace gini --column counts ctrl_counts.csv  # Expecting a low gini
outerspace gini --column exp_counts.csv counts # Expecting a high gini

# Find significant protospacers
outerspace diff --column counts ctrl_counts.csv exp_counts.csv results.csv 
```

## Barcoded Viruses for Latency Studies

Barcoded viruses, such as SIVmac293m2, are powerful tools for studying viral latency and reservoir dynamics in animal models.
These viruses contain unique molecular barcodes integrated into their genome, allowing researchers to track individual viral lineages.

Here's how they work:

1. A pool of viruses, each containing a unique barcode sequence, is used to infect the animal model
2. During infection, each virus integrates its genome (including the barcode) into host cells
3. Some infected cells become latently infected, harboring dormant virus
4. When sampling tissues, the barcodes can be sequenced to:
   - Identify which viral variants established latent infection
   - Track the clonal expansion of infected cells
   - Map the anatomical distribution of viral reservoirs
   - Monitor changes in the reservoir composition over time

### Using OUTERSPACE for Barcode Analysis

OUTERSPACE is well-suited for analyzing barcode sequencing data from these experiments:

1. **Sequence Extraction**: The `find` command can extract viral barcodes from sequencing reads using configurable patterns that match the barcode context.

2. **Error Correction**: The `collapse` command corrects sequencing errors in barcodes, ensuring accurate lineage tracking. This is crucial because even single nucleotide errors could artificially inflate diversity estimates.

3. **Quantification**: The `count` command determines the frequency of each viral barcode in different samples, revealing the relative abundance of viral variants.

4. **Distribution Analysis**: The `gini` command calculates Gini coefficients to measure the inequality of barcode distributions, helping identify:
   - Bottleneck effects during transmission
   - Clonal expansion of infected cells
   - Changes in reservoir diversity over time

5. **Comparative Analysis**: The `diff` command can compare barcode frequencies between:
   - Different anatomical sites
   - Pre- and post-treatment timepoints
   - Different experimental conditions

This integrated analysis pipeline helps researchers understand the dynamics of viral reservoirs, evaluate therapeutic interventions, and track viral dissemination throughout the host.

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
