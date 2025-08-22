"""Wrapper for outerspace merge command"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.2"

from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_files = snakemake.input.csv  # Only get the CSV files, not the TOML
output_file = snakemake.output[0]
toml_file = snakemake.input.get("toml", None)


# Construct command line arguments
args = [
    'merge',
    '--output-file', output_file,
    '-c', toml_file,
]

# Add input files (only CSV files)
args.extend(input_files)

# Run the merge command
cli = Cli(args)
cli.run() 