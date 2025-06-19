# Gini Wrapper

This wrapper provides a Snakemake interface for the outerspace gini command, which calculates the Gini coefficient from counts in a CSV column.

## Input/Output

### Input
- `input`: CSV file containing the data to analyze

### Output
- `output`: Text file containing the calculated Gini coefficient

### Parameters
- `column`: Column name to calculate Gini coefficient from (default: "counts")
- `count_column`: Optional column containing pre-counted values
- `scale`: Optional scale factor for normalized values (e.g., if normalized to mean=1)
- `sep`: CSV separator (default: ",")
- `allowed_list`: Optional text file containing allowed values, one per line

## Example

```python
rule gini:
    input:
        "path/to/input.csv"
    output:
        "path/to/gini.txt"
    params:
        column = "counts",
        count_column = "frequency",  # optional
        scale = 1000,  # optional
        allowed_list = "path/to/allowed.txt"  # optional
    wrapper:
        "wrappers/gini"
```

## Notes

- The input CSV must contain the specified column
- If count_column is provided, values will be weighted by their counts
- If scale is provided, counts will be multiplied by this factor
- If allowed_list is provided, only values in the list will be considered
- The output file will contain a single number representing the Gini coefficient
- The Gini coefficient ranges from 0 (perfect equality) to 1 (perfect inequality) 