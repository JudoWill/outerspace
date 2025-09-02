# Stats Wrapper

This wrapper provides a Snakemake interface for the outerspace stats command, which calculates all single-library statistics from counts in CSV files.

## Input/Output

### Input
- `input`: CSV file(s) containing UMI counts to analyze (can be single file or list of files)

### Output
- `output`: CSV file containing calculated statistics for all input files

### Parameters
- `key_column`: Column name containing keys/UMIs (default: "protospacer")
- `count_column`: Column name containing pre-counted values (optional)
- `scale`: Scale factor for normalized values (optional)
- `sep`: CSV separator (default: ",")
- `allowed_list`: Text file containing allowed values, one per line (optional)

## Example

```python
# Single file
rule stats_single:
    input:
        "counts.csv"
    output:
        "statistics.csv"
    params:
        key_column = "protospacer",
        count_column = "umi3_umi5_corrected_count",
        scale = 1.0,  # optional
        allowed_list = "allowed_values.txt"  # optional
    wrapper:
        "wrappers/stats"

# Multiple files
rule stats_multiple:
    input:
        ["sample1_counts.csv", "sample2_counts.csv", "sample3_counts.csv"]
    output:
        "all_statistics.csv"
    params:
        key_column = "protospacer",
        count_column = "umi3_umi5_corrected_count"
    wrapper:
        "wrappers/stats"
```

## Notes

- The input CSV file(s) must contain the specified key column
- If count_column is not provided, assumes count=1 for each key
- Scale factor can be used for normalized values (e.g., if normalized to mean=1)
- Allowed list filters calculations to only include specified values
- Output includes all calculated statistics: Gini coefficient, Shannon diversity, Simpson diversity, UMI recovery rate, UMI efficiency rate, UMI error rate, and UMI redundancy
- When multiple files are provided, statistics for each file are included as separate rows in the output CSV
- The output CSV includes a 'filename' column to identify which file each row of statistics corresponds to 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.