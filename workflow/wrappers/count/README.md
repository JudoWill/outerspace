# Count Wrapper

This wrapper provides a Snakemake interface for the outerspace count command, which counts unique barcodes per key value in CSV files.

## Input/Output

### Input
- `input`: CSV file containing barcodes and keys to count

### Output
- `output`: CSV file containing barcode counts per key

### Parameters
- `barcode_column`: Column name containing barcodes (default: "umi3_umi5_corrected")
- `key_column`: Column name to group by (default: "protospacer")
- `sep`: CSV separator (default: ",")
- `detailed`: Include barcode lists in output (default: false)
- `downsample`: Randomly sample reads with probability between 0 and 1 (optional)
- `random_seed`: Random seed for downsampling (optional)
- `allowed_list`: Text file containing allowed keys, one per line (optional)
- `metrics`: Path to output YAML file for metrics (optional)

## Example

```python
rule count:
    input:
        "path/to/input.csv"
    output:
        "path/to/output.csv"
    params:
        barcode_column = "umi3_umi5_corrected",
        key_column = "protospacer",
        detailed = True,
        downsample = 0.5,  # optional
        random_seed = 42,  # optional
        allowed_list = "path/to/allowed.txt",  # optional
        metrics = "path/to/metrics.yaml"  # optional
    wrapper:
        "wrappers/count"
```

## Notes

- The input CSV must contain the specified barcode and key columns
- The output CSV will contain key and barcode count columns
- If detailed output is enabled, the output will also include a list of unique barcodes for each key
- If an allowed list is provided, only keys in the list will be counted
- If metrics file is specified, it will contain statistics about the counting process
- Downsampling can be used to randomly sample reads with a specified probability 