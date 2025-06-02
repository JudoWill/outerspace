"""Wrapper for outerspace count barcode counting"""

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
barcode_column = snakemake.params.get("barcode_column", "umi3_umi5_corrected")
key_column = snakemake.params.get("key_column", "protospacer")
sep = snakemake.params.get("sep", ",")
detailed = snakemake.params.get("detailed", False)
downsample = snakemake.params.get("downsample", None)
random_seed = snakemake.params.get("random_seed", None)
allowed_list = snakemake.params.get("allowed_list", None)
metrics_file = snakemake.params.get("metrics", None)

# Construct command line arguments
args = [
    'count',
    '--input-file', input_file,
    '--output-file', output_file,
    '--barcode-column', barcode_column,
    '--key-column', key_column,
    '--sep', sep
]

if detailed:
    args.append('--detailed')
if downsample is not None:
    args.extend(['--downsample', str(downsample)])
if random_seed is not None:
    args.extend(['--random-seed', str(random_seed)])
if allowed_list:
    args.extend(['--allowed-list', allowed_list])
if metrics_file:
    args.extend(['--metrics', metrics_file])

# Run the count command
cli = Cli(args)
cli.run() 