"""Wrapper for outerspace collapse barcode correction"""

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
input_file = snakemake.input[0]
toml_file = snakemake.input.get("toml", None)
output_file = snakemake.output[0]

# Get parameters with defaults
columns = snakemake.params.get("columns", "umi3,umi5")
mismatches = snakemake.params.get("mismatches", 2)
method = snakemake.params.get("method", "directional")
sep = snakemake.params.get("sep", ",")
metrics_file = snakemake.params.get("metrics", None)

# Construct command line arguments
args = [
    'collapse',
    '--input-file', input_file,
    '--output-file', output_file,
    '--columns', columns,
    '--mismatches', str(mismatches),
    '--method', method,
    '--sep', sep
]

if metrics_file:
    args.extend(['--metrics', metrics_file])

if toml_file:
    args.extend(['--config', toml_file])

# Run the collapse command
cli = Cli(args)
cli.run() 