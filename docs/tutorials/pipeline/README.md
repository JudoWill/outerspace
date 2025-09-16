# OUTERSPACE Pipeline Tutorial

This tutorial demonstrates how to use the OUTERSPACE pipeline command with Snakemake to automate the complete analysis workflow. We'll use the same SIV barcoding data from the previous tutorial to show how pipeline automation can streamline large-scale analyses.

## Overview

The OUTERSPACE pipeline uses Snakemake to orchestrate the complete analysis workflow:
1. **findseq**: Extract sequences from BAM/FASTQ files
2. **collapse**: Correct barcodes using UMI-tools clustering  
3. **count**: Quantify unique barcodes per sample
4. **merge**: Combine results across samples
5. **stats**: Generate comprehensive statistics

Using the pipeline approach offers several advantages:
- **Automated execution**: All steps run automatically with dependency tracking
- **Parallelization**: Multiple samples processed simultaneously
- **Resume capability**: Failed jobs can be restarted without losing progress
- **Scalability**: Easy integration with HPC clusters and job schedulers
- **Reproducibility**: Consistent execution across different environments

## Prerequisites

1. Navigate to this tutorial directory:
```bash
cd docs/tutorials/pipeline
```

2. Copy SIV data for pipeline demonstration:
```bash
# Create data directory and copy BAM files
mkdir -p data
cp ../siv-barcoding/data/*.bam data/
cp ../siv-barcoding/sivbarcode.toml .
```

## Configuration Files

The pipeline requires two configuration files:

### 1. TOML Configuration (sivbarcode.toml)
This defines the patterns and command parameters (same as individual command tutorials):

```toml
# Global patterns for viral barcode extraction
[[patterns]]
name = "viral_barcode"
reg_expr = "(cccctccaggactagcataa){i<=1,d<=1,s<=2}(?P<viral_barcode>.{34})(atggaagaaagacctccaga){i<=1,d<=1,s<=2}"
read = "R1"
orientation = "both"
multiple = "first"

[findseq]
use_all_patterns = true
matches_only = true
threads = 8
skip_unmapped = true

[collapse]
columns = 'viral_barcode'
mismatches = 2

[count]
barcode_column = 'viral_barcode_corrected'
key_column = 'viral_barcode_corrected'

[merge]
key_column = 'viral_barcode_corrected'
count_column = 'viral_barcode_corrected_count'

[stats]
key_column = 'viral_barcode_corrected'
count_column = 'viral_barcode_corrected_count'
```

### 2. Snakemake Configuration (snakemake_config.yaml)
This defines the input data and pipeline parameters:

```yaml
# OUTERSPACE Pipeline Configuration for SIV Barcoding Analysis

# TOML configuration file containing patterns and command parameters
toml: "sivbarcode.toml"

# Method 1: Directory-based input (recommended for this tutorial)
direc: "data"
paired: false  # BAM files are single-end

# Method 2: Sample sheet (alternative approach)
# samples: "samplesheet.csv"

# Output sample names (optional - will use filenames if not specified)
sample_names:
  - "Early_2023-04-06"
  - "Later_2023-07-11"

# Pipeline execution parameters
cores: 4
executor: "local"  # Options: local, slurm, cluster

# Resource specifications (useful for cluster execution)
resources:
  findseq:
    threads: 2
    mem_mb: 4000
    time: "02:00:00"
  collapse:
    threads: 1
    mem_mb: 2000
    time: "01:00:00"
  count:
    threads: 1
    mem_mb: 1000
    time: "00:30:00"
  merge:
    threads: 1
    mem_mb: 1000
    time: "00:15:00"
  stats:
    threads: 1
    mem_mb: 1000
    time: "00:15:00"
```

## Basic Pipeline Execution

### 1. Dry Run (Testing)
Always start with a dry run to verify the workflow:

```bash
# Test the pipeline without executing
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--dry-run"
```

This shows:
- Which rules will be executed
- Input/output file dependencies
- Estimated resource requirements
- Any configuration errors

### 2. Local Execution
Run the pipeline on your local machine:

```bash
# Execute with 4 cores
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4"

# Execute with verbose output
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --verbose"

# Execute with progress tracking
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --progress"
```

### 3. Partial Execution
Run only specific parts of the pipeline:

```bash
# Only run findseq step
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 findseq_all"

# Run up to count step
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 count_all"

# Generate only statistics
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 stats_all"
```

## Advanced Snakemake Options

### Resource Management

```bash
# Limit memory usage per job
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --resources mem_mb=8000"

# Set maximum number of jobs
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --jobs 2"

# Use specific number of threads per rule
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --set-threads findseq=8"
```

### Debugging and Monitoring

```bash
# Print shell commands being executed
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --printshellcmds"

# Keep temporary files for debugging
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --notemp"

# Print detailed rule information
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --summary"

# Show rule graph
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --rulegraph | dot -Tpng > rulegraph.png"
```

### Restarting and Cleanup

```bash
# Force re-execution of all rules
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --forceall"

# Re-run specific samples
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --forcerun findseq"

# Clean up output files
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--delete-all-output"
```

## HPC Cluster Execution

### SLURM Integration

For large-scale analyses on HPC clusters with SLURM:

```bash
# Submit to SLURM cluster
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--executor slurm --jobs 100"

# With specific SLURM parameters
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--executor slurm --jobs 100 --default-resources slurm_account=myaccount slurm_partition=compute"

# Dry run for cluster submission
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--executor slurm --jobs 100 --dry-run"
```

### PBS/Torque Integration

```bash
# Submit to PBS cluster
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--executor cluster-generic --cluster 'qsub -V -cwd -pe smp {threads}' --jobs 50"
```

### Cluster Configuration Example

Create a cluster configuration file `cluster_config.yaml`:

```yaml
__default__:
  time: "04:00:00"
  mem: "4G"
  partition: "compute"
  account: "myproject"

findseq:
  time: "08:00:00"
  mem: "8G"
  threads: 8

collapse:
  time: "02:00:00"
  mem: "4G"
  threads: 2

count:
  time: "01:00:00"
  mem: "2G"
  threads: 1
```

Then run with:
```bash
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--executor slurm --cluster-config cluster_config.yaml --jobs 100"
```

## Sample Sheet Alternative

Instead of using directory-based input, you can use a sample sheet for more control:

Create `samplesheet.csv`:
```csv
sample_name,reads,toml
Early_Timepoint,data/MJ21-20230406.sifted.bam,sivbarcode.toml
Later_Timepoint,data/MJ21-20230711.sifted.bam,sivbarcode.toml
```

Update `snakemake_config.yaml`:
```yaml
toml: "sivbarcode.toml"
samples: "samplesheet.csv"
```

This approach allows:
- Custom sample naming
- Different TOML configs per sample
- Paired-end file specification
- Metadata inclusion

## Pipeline Outputs

The pipeline creates a structured output directory:

```
├── findseq/
│   ├── Early_Timepoint.csv
│   └── Later_Timepoint.csv
├── collapse/
│   ├── Early_Timepoint.csv
│   └── Later_Timepoint.csv
├── count/
│   ├── Early_Timepoint.csv
│   └── Later_Timepoint.csv
├── merge/
│   └── merged.csv
└── stats/
    └── stats.csv
```

## Monitoring and Logs

Snakemake provides comprehensive logging:

```bash
# View recent logs
ls -la .snakemake/log/

# Monitor job status
snakemake --summary

# Check failed jobs
snakemake --detailed-summary | grep FAILED
```

## Best Practices

### 1. Development Workflow
```bash
# Start with dry run
outerspace pipeline sivbarcode.toml snakemake_config.yaml --snakemake-args="--dry-run"

# Test with single core
outerspace pipeline sivbarcode.toml snakemake_config.yaml --snakemake-args="--cores 1"

# Scale up gradually
outerspace pipeline sivbarcode.toml snakemake_config.yaml --snakemake-args="--cores 4"

# Deploy to cluster
outerspace pipeline sivbarcode.toml snakemake_config.yaml --snakemake-args="--executor slurm --jobs 100"
```

### 2. Resource Optimization
- Start with conservative resource estimates
- Monitor actual usage with `--benchmark-repeats`
- Adjust based on profiling results
- Use cluster-specific optimizations

### 3. Error Handling
- Always use `--dry-run` first
- Keep intermediate files during development (`--notemp`)
- Use `--printshellcmds` for debugging
- Monitor cluster queue status

### 4. Reproducibility
- Version control all configuration files
- Document resource requirements
- Use container execution when possible
- Archive complete workflows

## Common Snakemake Arguments Reference

| Argument | Purpose | Example |
|----------|---------|---------|
| `--dry-run` or `-n` | Test without execution | `--snakemake-args="--dry-run"` |
| `--cores` or `-c` | Set CPU cores | `--snakemake-args="--cores 8"` |
| `--jobs` or `-j` | Max parallel jobs | `--snakemake-args="--jobs 100"` |
| `--executor` | Execution backend | `--snakemake-args="--executor slurm"` |
| `--forceall` | Re-run everything | `--snakemake-args="--forceall"` |
| `--forcerun` | Re-run specific rule | `--snakemake-args="--forcerun findseq"` |
| `--printshellcmds` | Show commands | `--snakemake-args="--printshellcmds"` |
| `--verbose` | Detailed output | `--snakemake-args="--verbose"` |
| `--summary` | Job summary | `--snakemake-args="--summary"` |
| `--delete-all-output` | Clean outputs | `--snakemake-args="--delete-all-output"` |

## Troubleshooting

### Common Issues

1. **Configuration errors**: Use `--dry-run` to catch early
2. **Resource constraints**: Monitor with `--summary` and adjust
3. **Failed jobs**: Check `.snakemake/log/` for details
4. **Dependency issues**: Verify input file paths
5. **Cluster problems**: Test locally first

### Performance Tips

1. **Parallel execution**: Use appropriate `--cores` and `--jobs`
2. **Resource allocation**: Match cluster specifications
3. **Input/output optimization**: Use fast storage for temporary files
4. **Memory management**: Monitor and adjust resource requests

## Example: Complete Pipeline Run

Here's a complete example with monitoring:

```bash
# 1. Verify configuration
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--dry-run --verbose"

# 2. Execute with monitoring
outerspace pipeline sivbarcode.toml snakemake_config.yaml \
    --snakemake-args="--cores 4 --progress --printshellcmds"

# 3. Check results
ls -la */
head merge/merged.csv
head stats/stats.csv
```

This tutorial demonstrates how the pipeline approach can dramatically simplify large-scale OUTERSPACE analyses while providing the flexibility and power of Snakemake workflow management.

Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
