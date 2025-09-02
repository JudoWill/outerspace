# FindSeq Wrapper

This wrapper provides a Snakemake interface for the outerspace findseq command, which extracts sequences from FASTQ files based on configuration patterns.

## Input/Output

### Input
- `reads`: Either a single FASTQ file or a list of two FASTQ files for paired-end sequencing
- `toml`: Configuration file containing search patterns

### Output
- `output`: CSV file containing the extracted sequences

## Example

```python
rule findseq:
    input:
        reads = "path/to/reads.fastq.gz",  # or ["path/to/read1.fastq.gz", "path/to/read2.fastq.gz"]
        toml = "path/to/config.toml"
    output:
        "path/to/output.csv"
    wrapper:
        "wrappers/findseq"
```

## Notes

- For paired-end sequencing, provide both read files as a list in the `reads` input
- The configuration file should be in TOML format and contain the search patterns
- The output CSV will contain columns: read_id, UMI_5prime, protospacer, downstreamof_protospacer, UMI_3prime 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.