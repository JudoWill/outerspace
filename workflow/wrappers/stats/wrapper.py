"""Wrapper for outerspace stats command"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.1"

import sys
from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_files = snakemake.input
output_file = snakemake.output[0]
toml_file = snakemake.input.get("toml", None)

# Get parameters with defaults
key_column = snakemake.params.get("key_column", "protospacer")
count_column = snakemake.params.get("count_column", None)
scale = snakemake.params.get("scale", None)
sep = snakemake.params.get("sep", ",")
allowed_list = snakemake.params.get("allowed_list", None)

# Construct command line arguments
args = [
    'stats',
    '--key-column', key_column,
    '--sep', sep
]

# Add input files (can be single file or list of files)
if isinstance(input_files, str):
    args.append(input_files)
else:
    args.extend(input_files)

if toml_file:
    args.extend(['--config', toml_file])

if count_column:
    args.extend(['--count-column', count_column])
if scale is not None:
    args.extend(['--scale', str(scale)])
if allowed_list:
    args.extend(['--allowed-list', allowed_list])

# Redirect stdout to output file since stats command writes to stdout
original_stdout = sys.stdout
with open(output_file, 'w') as f:
    sys.stdout = f
    try:
        # Run the stats command
        cli = Cli(args)
        cli.run()
    finally:
        # Restore stdout
        sys.stdout = original_stdout 
        
        
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.