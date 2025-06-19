"""Wrapper for outerspace findseq sequence extraction"""

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
config_file = snakemake.input.toml
read1_file = snakemake.input.reads[0] if isinstance(snakemake.input.reads, list) else snakemake.input.reads
read2_file = snakemake.input.reads[1] if isinstance(snakemake.input.reads, list) else None
output_file = snakemake.output[0]

# Construct command line arguments
args = ['findseq', config_file, '-1', read1_file]
if read2_file:
    args.extend(['-2', read2_file])
args.extend(['-o', output_file])

# Run the findseq command
cli = Cli(args)
cli.run() 