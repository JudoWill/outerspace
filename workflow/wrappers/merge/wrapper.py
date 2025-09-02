"""Wrapper for outerspace merge command"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.1"

from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_files = snakemake.input.csv  # Only get the CSV files, not the TOML
output_file = snakemake.output[0]
toml_file = snakemake.input.get("toml", None)

# Get parameters with defaults
key_column = snakemake.params.get("key_column", "protospacer")
count_column = snakemake.params.get("count_column", None)
sample_names = snakemake.params.get("sample_names", None)
sep = snakemake.params.get("sep", ",")
format_type = snakemake.params.get("format", "wide")
mismatches = snakemake.params.get("mismatches", 0)
method = snakemake.params.get("method", "directional")
metrics_file = snakemake.params.get("metrics", None)

# Construct command line arguments
args = [
    'merge',
    '--output-file', output_file,
    '--key-column', key_column,
    '--sep', sep,
    '--format', format_type,
    '--mismatches', str(mismatches),
    '--method', method
]

# Add input files (only CSV files)
args.extend(input_files)

if toml_file:
    args.extend(['--config', toml_file])

if count_column:
    args.extend(['--count-column', count_column])
if sample_names:
    args.extend(['--sample-names'] + sample_names)
if metrics_file:
    args.extend(['--metrics', metrics_file])

# Run the merge command
cli = Cli(args)
cli.run() 

# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.