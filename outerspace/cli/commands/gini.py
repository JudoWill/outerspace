"""Gini coefficient calculation command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import csv
import sys
from typing import List, Optional
from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.stats import GiniCoefficient

class GiniCommand(BaseCommand):
    """Command for calculating Gini coefficient from counts in a CSV column"""
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('gini',
            help='Calculate Gini coefficient from counts in a CSV column')
        parser.add_argument('input_file',
            help='Input CSV file')
        parser.add_argument('--column', required=True,
            help='Column to calculate Gini coefficient from')
        parser.add_argument('--count-column',
            help='Column containing pre-counted values')
        parser.add_argument('--scale', type=float,
            help='Scale factor for normalized values (e.g., if normalized to mean=1)')
        parser.add_argument('--sep', default=',',
            help='CSV separator (default: ,)')
        parser.add_argument('--allowed-list',
            help='Text file containing allowed values (one per line)')
        parser.add_argument('--config',
            help='TOML configuration file containing command settings')
        return parser

    def _read_allowed_list(self, filepath: str) -> List[str]:
        """Read allowed values from a text file"""
        with open(filepath, 'r') as f:
            return [line.strip() for line in f if line.strip()]

    def _calculate_gini(self, input_file: str, column: str, count_column: Optional[str] = None, 
                       scale: Optional[float] = None, sep: str = ",", 
                       allowed_list: Optional[List[str]] = None) -> Optional[float]:
        """Calculate Gini coefficient for counts in a CSV column
        
        Args:
            input_file: Path to input CSV file
            column: Column name to calculate Gini coefficient from
            count_column: Optional column containing pre-counted values
            scale: Optional scale factor for normalized values
            sep: CSV separator
            allowed_list: Optional list of allowed values
            
        Returns:
            Gini coefficient or None if no data available
        """
        # Create UMI object from CSV
        umi = UMI.from_csv(
            filepath=input_file,
            column=column,
            mismatches=0,  # No correction needed for Gini calculation
            method="adjacency",  # Method doesn't matter since mismatches=0
            sep=sep,
            correct=False,  # No correction needed for Gini calculation
            count_column=count_column,
            scale=scale
        )
        
        # Calculate Gini coefficient using the stats module
        result = GiniCoefficient.calculate(umi, allowed_list=allowed_list)
        gini = result['gini_coefficient']
        
        # Get some statistics for reporting
        counts = umi.corrected_counts
        total_unique = len(counts)
        total_counts = sum(counts.values())
        
        # Print statistics
        print(f"\nStatistics:", file=sys.stderr)
        print(f"  Total unique values: {total_unique}", file=sys.stderr)
        print(f"  Total counts: {total_counts}", file=sys.stderr)
        if scale:
            print(f"  Scale factor applied: {scale}", file=sys.stderr)
        print(f"  Gini coefficient: {gini:.4f}" if gini is not None else "  Gini coefficient: None", file=sys.stderr)
        
        if allowed_list:
            missing = set(allowed_list) - {k.decode('ascii') for k in counts.keys()}
            print(f"  Missing values: {len(missing)}/{len(allowed_list)}", file=sys.stderr)
            if missing and len(missing) <= 10:
                print("  Missing values list:", file=sys.stderr)
                for value in sorted(missing):
                    print(f"    {value}", file=sys.stderr)
        
        return gini

    def run(self):
        """Run the gini command"""
        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)
        
        # Merge config and args with defaults
        defaults = {
            'sep': ',',
            'scale': None
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.input_file:
            raise ValueError("Please provide an input file")
        if not self.args.column:
            raise ValueError("Please provide a column to calculate Gini coefficient for")

        # Validate scale if provided
        if self.args.scale is not None and self.args.scale <= 0:
            raise ValueError("Scale value must be positive")

        try:
            # Read allowed list if provided
            allowed_list = None
            if self.args.allowed_list:
                allowed_list = self._read_allowed_list(self.args.allowed_list)
                print(f"Loaded {len(allowed_list)} allowed values from {self.args.allowed_list}", file=sys.stderr)
            
            # Calculate Gini coefficient
            gini = self._calculate_gini(
                self.args.input_file,
                self.args.column,
                count_column=self.args.count_column,
                scale=self.args.scale,
                sep=self.args.sep,
                allowed_list=allowed_list
            )
            
            # Print result to stdout (for scripting)
            if gini is not None:
                print(f"{gini:.6f}")
            else:
                print("NA")
                sys.exit(1)
                
        except Exception as e:
            raise ValueError(f"Error calculating Gini coefficient: {e}") 