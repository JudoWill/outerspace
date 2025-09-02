# Merge Wrapper

This wrapper provides a Snakemake interface for the outerspace merge command, which merges multiple UMI count files into a single file.

## Input/Output

### Input
- `input`: List of CSV files containing UMI counts to merge

### Output
- `output`: CSV file containing merged UMI counts

### Parameters
- `key_column`: Column name containing UMIs (default: "protospacer")
- `count_column`: Column name containing counts (optional, assumes count=1 if not provided)
- `sample_names`: List of sample names matching input files (optional)
- `sep`: CSV separator (default: ",")
- `format`: Output format - "wide" (samples as columns) or "long" (sample,umi,count columns) (default: "wide")
- `mismatches`: Number of mismatches allowed for clustering (default: 0)
- `method`: Clustering method - "cluster", "adjacency", or "directional" (default: "directional")
- `metrics`: Path to output YAML file for metrics (optional)

## Example

```python
rule merge:
    input:
        ["sample1.csv", "sample2.csv", "sample3.csv"]
    output:
        "merged_counts.csv"
    params:
        key_column = "protospacer",
        count_column = "umi3_umi5_corrected_count",
        sample_names = ["sample1", "sample2", "sample3"],
        format = "wide",
        mismatches = 2,
        method = "directional",
        metrics = "merge_metrics.yaml"  # optional
    wrapper:
        "wrappers/merge"
```

## Notes

- The input CSV files must contain the specified key column
- If count_column is not provided, assumes count=1 for each UMI
- Sample names must match the number of input files if provided
- Wide format creates columns for each sample, long format creates sample,umi,count columns
- If mismatches > 0, UMIs will be clustered before merging
- Metrics file will contain statistics about the merging and clustering process 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.