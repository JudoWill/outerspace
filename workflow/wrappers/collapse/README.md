# Collapse Wrapper

This wrapper provides a Snakemake interface for the outerspace collapse command, which corrects barcodes in CSV files using UMI-tools clustering.

## Input/Output

### Input
- `input`: CSV file containing barcodes to correct

### Output
- `output`: CSV file containing corrected barcodes

### Parameters
- `columns`: Comma-separated list of column names containing barcodes to correct (default: "umi3,umi5")
- `mismatches`: Number of mismatches allowed for clustering (default: 2)
- `method`: Clustering method to use (default: "directional", options: "cluster", "adjacency", "directional")
- `sep`: CSV separator (default: ",")
- `metrics`: Optional path to output YAML file for metrics

## Example

```python
rule collapse:
    input:
        "path/to/input.csv"
    output:
        "path/to/output.csv"
    params:
        columns = "umi3,umi5",
        mismatches = 2,
        method = "directional",
        metrics = "path/to/metrics.yaml"  # optional
    wrapper:
        "wrappers/collapse"
```

## Notes

- The input CSV must contain the specified columns
- The output CSV will contain all original columns plus a new column with the corrected barcodes
- The corrected column name will be the concatenation of the input columns with "_corrected" suffix
- If metrics file is specified, it will contain statistics about the correction process 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.