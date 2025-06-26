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

## Global Patterns

OUTERSPACE now uses a global patterns system where patterns are defined once and can be reused across multiple commands:

```toml
# Global patterns that can be used across multiple commands
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
```

### Pattern Configuration

Each pattern requires the following fields:

- `name`: Unique identifier for the pattern (required for global patterns)
- `reg_expr`: Regular expression with named capture groups
- `read`: Which read to search ("R1", "R2", or "both")
- `orientation`: Search orientation ("forward", "reverse-complement", or "both")
- `multiple`: How to handle multiple matches ("first", "last", or "all")
- `left_flank`: Number of letters to include left of the match. Defaults to 0.
- `right_flank`: Number of letters to include right of the match. Defaults to 0.

## Command Sections

Each command has its own section in the configuration file. The available sections are:

### [findseq]
Configuration for sequence extraction:
```toml
[findseq]
# Reference patterns by name
pattern_names = ["UMI_5prime", "protospacer", "UMI_3prime"]

# Or use all global patterns
# use_all_patterns = true
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

[collapse]
input_file = "results/matches.csv"
output_file = "results/corrected.csv"
columns = "UMI_5prime,UMI_3prime"
mismatches = 2
method = "directional"

[count]
input_file = "results/corrected.csv"
output_file = "results/counts.csv"
barcode_column = "UMI_5prime_UMI_3prime_corrected"
key_column = "protospacer"
detailed = true

[gini]
input_file = "results/counts.csv"
column = "UMI_5prime_UMI_3prime_corrected_count"
scale = 1.0

[pipeline]
config_file = "config.toml"
snakemake_config = "snakemake.yaml"
snakemake_args = "--cores 4"
```

## Best Practices

1. Use descriptive pattern names
2. Group related patterns in the global patterns section
3. Document non-default values with comments
4. Use relative paths when possible
5. Keep sensitive information out of configuration files
6. Version control your configuration files
7. Use separate configuration files for different projects or experiments

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
