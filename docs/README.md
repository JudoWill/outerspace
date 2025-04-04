# OuterSpace

Outerspace is a collection of tools for analyzing pooled CRISPR screens, viral barcode population studies, and any other application that requires the extraction and couting of variable regions in pooled amplicons.
It contains tools to extract regions of interest, correct sequencing error, assess diversity, and compare between samples.

## Contents
- [Installation](install.md)
- [Quick Start](#quick-start)
- [Biologic Applications](usage.md)
- [Basic Commands](commands.md)
- [Configuration](config.md)

## Quick Start

### Install

```bash
pip install outerspace
```

### Design your extraction strategy

`outerspace` uses the `regex` library to extract relevant features from a DNA sequence.
This allows an simple, expressive, and modular strategy for extracting of regions of interest while tolerating mismatches.
It supports both short, paired end reads and log long reads.
See the [walkthrough](regex_explainer.md) for a detailed discussion on how to design your extraction strategy.

### Create your config file

If you are going to repeating similar experiments often, `outerspace` allows you to encapsulate that information in a `toml` file accepted by all commands.
This ensures repeatability between analyses and can drastically simplify analyses.
It also facilitates reproducible science as the config can be stored, shared, and tracked.

See the [walkthrough](config.md) for a more detailed discussion on creating your config file.

### Commands

Assuming you are performing a CRISPR screen where you have used Illumina paired-end sequencing and you have sequenced your library before (_ctrl_) and after (_exp_) selection as described in the [walkthrough](regex_explainer.md).

```bash
# Extract sequence information
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


Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
