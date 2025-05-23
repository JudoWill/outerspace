# OuterSpace

Outerspace is a collection of tools for analyzing pooled CRISPR screens, viral barcode population studies, and any other application that requires the extraction and couting of variable regions in pooled amplicons.
It contains tools to extract regions of interest, correct sequencing error, assess diversity, and compare between samples.

## Contents
- [Installation](install.md)
- [Quick Start](#quick-start)
- [Biologic Applications](usage.md)
- [Basic Commands](commands.md)
- [Configuration](config.md)
- [Detailed Walkthrough](walkthrough.md)

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

### Process Your Data

For most analyses, you can use the `pipeline` command to process all your data in one step:

```bash
# Create output directory
mkdir -p results

# Run the pipeline
outerspace pipeline config.toml \
    --input-dir fastq_files \
    --output-dir results \
    --barcode-columns UMI_5prime,UMI_3prime \
    --key-column protospacer \
    --mismatches 2 \
    --method directional \
    --metrics
```

This will:
1. Process all FASTQ files in your input directory
2. Extract sequences using your config patterns
3. Correct barcodes using UMI-tools clustering
4. Count unique barcodes per protospacer
5. Generate metrics files for quality control

For more detailed instructions, including how to run individual commands and perform additional analyses, see the [detailed walkthrough](walkthrough.md).

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
