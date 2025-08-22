"""Wrapper for outerspace stats command"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.2"

import sys
from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_files = snakemake.input
output_file = snakemake.output[0]
toml_file = snakemake.input.get("toml", None)

# Construct command line arguments
args = [
    'stats',
    '-c', toml_file,
]

# Add input files (can be single file or list of files)
if isinstance(input_files, str):
    args.append(input_files)
else:
    args.extend(input_files)

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