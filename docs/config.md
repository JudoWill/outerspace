# Configuration

OUTERSPACE uses TOML (Tom's Obvious, Minimal Language) configuration files to set parameters for its commands. This document describes how to use the configuration files effectively.

## TOML Syntax

TOML is a simple configuration file format that's easy to read and write. Here are the key syntax elements used in OUTERSPACE:

```toml
# Comments start with #
[section]  # Section headers in square brackets
key = "value"  # String values in quotes
number = 42  # Numbers without quotes
boolean = true  # Boolean values
array = ["item1", "item2"]  # Arrays in square brackets
```

## Command Sections

Each command has its own section in the configuration file. The available sections are:

### [findseq]
Configuration for sequence extraction:
```toml
[findseq]
config = "path/to/config.toml"  # Required
read1_filename = "path/to/read1.fastq.gz"  # Optional
read2_filename = "path/to/read2.fastq.gz"  # Optional
output_filename = "path/to/output.csv"  # Optional
fastqfiles = ["dir1", "dir2"]  # Optional
outdir = "path/to/output"  # Optional
read_regxlist = "pattern1,pattern2"  # Optional
read1_regxlist = "pattern1,pattern2"  # Optional
read2_regxlist = "pattern1,pattern2"  # Optional
```

### [collapse]
Configuration for barcode correction:
```toml
[collapse]
input_file = "path/to/input.csv"  # Optional
input_dir = "path/to/input"  # Optional
output_file = "path/to/output.csv"  # Optional
output_dir = "path/to/output"  # Optional
columns = "col1,col2"  # Required
mismatches = 2  # Optional, default: 2
sep = ","  # Optional, default: ","
row_limit = 1000  # Optional
method = "directional"  # Optional, default: "directional"
metrics = "path/to/metrics.yaml"  # Optional
```

### [count]
Configuration for barcode counting:
```toml
[count]
input_file = "path/to/input.csv"  # Optional
input_dir = "path/to/input"  # Optional
output_file = "path/to/output.csv"  # Optional
output_dir = "path/to/output"  # Optional
barcode_column = "barcode"  # Required
key_column = "key"  # Required
sep = ","  # Optional, default: ","
row_limit = 1000  # Optional
allowed_list = "path/to/allowed.txt"  # Optional
detailed = false  # Optional, default: false
downsample = 0.5  # Optional
random_seed = 42  # Optional
metrics = "path/to/metrics.yaml"  # Optional
```

### [gini]
Configuration for Gini coefficient calculation:
```toml
[gini]
input_file = "path/to/input.csv"  # Required
column = "counts"  # Required
count_column = "precounts"  # Optional
scale = 1.0  # Optional
sep = ","  # Optional, default: ","
allowed_list = "path/to/allowed.txt"  # Optional
```

### [pipeline]
Configuration for the complete pipeline:
```toml
[pipeline]
config_file = "path/to/config.toml"  # Required
snakemake_config = "path/to/snakemake.yaml"  # Required
snakemake_args = "--dry-run --cores 4"  # Optional
```

## Using Configuration Files

1. Create a TOML file with the desired sections and parameters
2. Pass the configuration file to the command using the `--config` option:
   ```bash
   outerspace collapse --config config.toml
   ```

3. Command-line arguments override configuration file values
4. Required parameters must be specified either in the config file or on the command line

## Example Configuration

Here's a complete example configuration file:

```toml
[findseq]
config = "search_patterns.toml"
read1_filename = "data/reads_R1.fastq.gz"
read2_filename = "data/reads_R2.fastq.gz"
output_filename = "results/matches.csv"

[collapse]
input_file = "results/matches.csv"
output_file = "results/corrected.csv"
columns = "UMI_5prime,UMI_3prime"
mismatches = 2
method = "directional"

[count]
input_file = "results/corrected.csv"
output_file = "results/counts.csv"
barcode_column = "UMI_5prime"
key_column = "protospacer"
detailed = true

[gini]
input_file = "results/counts.csv"
column = "count"
scale = 1.0

[pipeline]
config_file = "config.toml"
snakemake_config = "snakemake.yaml"
snakemake_args = "--cores 4"
```

## Best Practices

1. Use descriptive file paths
2. Group related parameters in the same section
3. Document non-default values with comments
4. Use relative paths when possible
5. Keep sensitive information out of configuration files
6. Version control your configuration files
7. Use separate configuration files for different projects or experiments

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
