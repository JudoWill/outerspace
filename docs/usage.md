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

1. **Sequence Extraction** (`findseq`): Extracts gRNA sequences and associated barcodes from FASTQ files using configurable search patterns.

2. **Barcode Correction** (`collapse`): Corrects sequencing errors in barcodes using UMI-tools algorithms to improve accuracy.

3. **Counting** (`count`): Quantifies the frequency of each gRNA-barcode combination.

4. **Distribution Analysis** (`gini`): Calculates Gini coefficients to measure the inequality of barcode distributions. This is useful for diagnostic purposes.

5. **Visualization** (`visualize`): Creates plots to help interpret the results.

This integrated toolset helps researchers accurately measure changes in gRNA abundance between conditions, identify hits from their screens, and ensure data quality through multiple analysis steps.

### Usage Story

```bash
# Extract sequences
outerspace findseq config.toml -1 ctrl_r1.fastq.gz -2 ctrl_r2.fastq.gz -o ctrl_output.csv
outerspace findseq config.toml -1 exp_r1.fastq.gz -2 exp_r2.fastq.gz -o exp_output.csv

# Single file processing
outerspace findseq config.toml -1 single_reads.fasta -o single_output.csv

# SAM/BAM with region specification
outerspace findseq config.toml -1 aligned_reads.bam --region "chr1:1-1000" -o region_output.csv

# Correct barcodes
outerspace collapse --input-file ctrl_output.csv --output-file ctrl_output_collapsed.csv \
    --columns umi3,umi5 --mismatches 2 --method directional
outerspace collapse --input-file exp_output.csv --output-file exp_output_collapsed.csv \
    --columns umi3,umi5 --mismatches 2 --method directional

# Count barcodes
outerspace count --input-file ctrl_output_collapsed.csv --output-file ctrl_counts.csv \
    --barcode-column umi3_umi5_corrected --key-column protospacer
outerspace count --input-file exp_output_collapsed.csv --output-file exp_counts.csv \
    --barcode-column umi3_umi5_corrected --key-column protospacer

# Calculate Gini coefficient
outerspace gini ctrl_counts.csv --column counts  # Expecting a low gini
outerspace gini exp_counts.csv --column counts  # Expecting a high gini

# Visualize results
outerspace visualize ctrl_counts.csv exp_counts.csv output_plots/ \
    --title-prefix "CRISPR Screen" --log-scale
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

1. **Sequence Extraction**: The `findseq` command can extract viral barcodes from sequencing reads using configurable patterns that match the barcode context.

2. **Error Correction**: The `collapse` command corrects sequencing errors in barcodes, ensuring accurate lineage tracking. This is crucial because even single nucleotide errors could artificially inflate diversity estimates.

3. **Quantification**: The `count` command determines the frequency of each viral barcode in different samples, revealing the relative abundance of viral variants.

4. **Distribution Analysis**: The `gini` command calculates Gini coefficients to measure the inequality of barcode distributions, helping identify:
   - Bottleneck effects during transmission
   - Clonal expansion of infected cells
   - Changes in reservoir diversity over time

5. **Visualization**: The `visualize` command creates plots to help interpret the distribution of viral barcodes across different samples.

### Example Pipeline

For analyzing viral barcode data, you can use the `pipeline` command to run all steps automatically:

```bash
outerspace pipeline config.toml snakemake_config.yaml \
    --snakemake-args="--cores 4"
```

This will:
1. Process all FASTQ files specified in the Snakemake configuration
2. Extract and correct barcodes using patterns in `config.toml`
3. Count unique barcodes per sample
4. Generate metrics and visualizations
5. Save all results in the output directory

Alternatively, you can run individual commands:

```bash
# Extract viral barcodes
outerspace findseq config.toml -1 viral_reads.fastq.gz -o viral_barcodes.csv

# Correct barcode errors
outerspace collapse --input-file viral_barcodes.csv --output-file corrected_barcodes.csv \
    --columns UMI_5prime,UMI_3prime --mismatches 2 --method directional

# Count barcodes per sample
outerspace count --input-file corrected_barcodes.csv --output-file barcode_counts.csv \
    --barcode-column UMI_5prime_UMI_3prime_corrected --key-column viral_barcode

# Analyze distribution
outerspace gini barcode_counts.csv --column UMI_5prime_UMI_3prime_corrected_count
```

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
