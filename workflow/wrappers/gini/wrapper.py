"""Wrapper for outerspace gini coefficient calculation"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd@example.com"
__license__ = "MIT"
__version__ = "0.0.1"

import os
from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Get parameters with defaults
column = snakemake.params.get("column", "counts")
count_column = snakemake.params.get("count_column", None)
scale = snakemake.params.get("scale", None)
sep = snakemake.params.get("sep", ",")
allowed_list = snakemake.params.get("allowed_list", None)

# Construct command line arguments
args = [
    'gini',
    input_file,
    '--column', column,
    '--sep', sep
]

if count_column:
    args.extend(['--count-column', count_column])
if scale is not None:
    args.extend(['--scale', str(scale)])
if allowed_list:
    args.extend(['--allowed-list', allowed_list])

# Run the gini command
cli = Cli(args)
cli.run()

# Write the output
with open(output_file, 'w') as f:
    f.write(str(cli.command._calculate_gini(
        input_file,
        column,
        count_column=count_column,
        scale=scale,
        sep=sep,
        allowed_list=allowed_list
    ))) 