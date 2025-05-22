"""Barcode counting command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from argparse import ArgumentParser
from sys import exit as sys_exit

from outerspace.cli.commands.base import BaseCommand

class CountCommand(BaseCommand):
    """Command for counting unique barcodes per key value in CSV files"""
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('count',
            help='Count unique barcodes per key value in CSV files')
        parser.add_argument('input_dir',
            help='Input directory containing CSV files')
        parser.add_argument('output_dir',
            help='Output directory for barcode counts')
        parser.add_argument('--barcode-column',
            help='Column containing barcodes')
        parser.add_argument('--key-column',
            help='Column to group by')
        parser.add_argument('--sep', default=',',
            help='CSV separator (default: ,)')
        parser.add_argument('--row-limit', type=int,
            help='Process only the first N rows (for testing)')
        parser.add_argument('--allowed-list',
            help='Text file containing allowed keys (one per line)')
        parser.add_argument('--detailed', action='store_true',
            help='Include barcode lists in output')
        parser.add_argument('--downsample', type=float,
            help='Randomly sample reads with probability between 0 and 1')
        parser.add_argument('--random-seed', type=int,
            help='Random seed for downsampling')
        parser.add_argument('--metrics',
            help='Output YAML file for metrics')
        return parser

    def run(self):
        """Run the count command"""
        # TODO: Implement barcode counting
        raise NotImplementedError("Barcode counting not yet implemented") 