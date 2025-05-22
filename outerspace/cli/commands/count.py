"""Barcode counting command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import csv
import glob
import os
import random
import sys
from collections import defaultdict
from typing import Dict, Set, Any
from tqdm import tqdm
import yaml

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI

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
        parser.add_argument('--barcode-column', required=True,
            help='Column containing barcodes')
        parser.add_argument('--key-column', required=True,
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

    def _read_allowed_keys(self, filepath: str) -> Set[str]:
        """Read allowed keys from a text file"""
        allowed_keys = set()
        with open(filepath, 'r') as f:
            for line in f:
                key = line.strip()
                if key:  # Skip empty lines
                    allowed_keys.add(key)
        return allowed_keys

    def _process_single_file(self, input_file: str, output_file: str, barcode_col: str, key_col: str,
                           sep: str, row_limit: int, allowed_keys: Set[str], detailed: bool,
                           downsample: float = None) -> Dict[str, Any]:
        """Process a single CSV file and return summary statistics"""
        # Create UMI object for this file
        umi = UMI(mismatches=0)
        
        # Read rows and collect barcodes per key
        barcodes_by_key = defaultdict(set)
        total_rows = 0
        rows_with_allowed_key = 0
        
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames
            
            # Verify columns exist
            missing_cols = [col for col in [barcode_col, key_col] if col not in headers]
            if missing_cols:
                raise ValueError(f"Columns not found in input file: {', '.join(missing_cols)}")
            
            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break
                    
                # Apply downsampling if specified
                if downsample is not None and random.random() > downsample:
                    continue
                    
                total_rows += 1
                key = str(row[key_col])
                barcode = str(row[barcode_col])
                
                # Skip if key is not in allowed list
                if allowed_keys:
                    if key in allowed_keys:
                        rows_with_allowed_key += 1
                        if key and barcode:  # Skip empty values
                            barcodes_by_key[key].add(barcode)
                            umi.consume(barcode)
                else:
                    if key and barcode:  # Skip empty values
                        barcodes_by_key[key].add(barcode)
                        umi.consume(barcode)
        
        # Calculate summary statistics
        total_keys = len(barcodes_by_key)
        total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())
        
        # Calculate Gini coefficient for barcodes (using UMI class)
        barcode_gini = umi.gini_coefficient()
        
        # Calculate Gini coefficient for keys using a new UMI instance
        key_umi = UMI(mismatches=0)  # No mismatches since keys are already corrected
        for key, barcodes in barcodes_by_key.items():
            key_umi.consume(key)
        key_umi.create_mapping()
        key_gini = key_umi.gini_coefficient(allowed_list=list(allowed_keys) if allowed_keys else None)
        
        stats = {
            'file_stats': {
                'total_rows_scanned': total_rows,
                'total_keys': total_keys,
                'total_barcodes': total_barcodes,
                'average_barcodes_per_key': total_barcodes / total_keys if total_keys > 0 else 0,
                'barcode_gini_coefficient': barcode_gini,
                'key_gini_coefficient': key_gini
            }
        }
        
        # Add allowed list statistics if allowed_keys is provided
        if allowed_keys:
            stats['file_stats']['rows_with_allowed_key'] = rows_with_allowed_key
            missing_keys = allowed_keys - set(barcodes_by_key.keys())
            stats['missing_keys'] = {
                'total_missing': len(missing_keys),
                'missing_keys_list': sorted(list(missing_keys))
            }
        
        # Write output
        self._write_counts(barcodes_by_key, output_file, sep, detailed)
        
        return stats

    def _write_counts(self, barcodes_by_key: Dict[str, Set[str]], filepath: str, sep: str, detailed: bool = False):
        """Write barcode counts per key to CSV file"""
        with open(filepath, 'w', newline='') as f:
            if detailed:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow(['key', 'unique_barcodes', 'barcode_count'])
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, ','.join(sorted(barcodes)), len(barcodes)])
            else:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow(['key', 'barcode_count'])
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, len(barcodes)])

    def _write_metrics(self, metrics: Dict[str, Any], filepath: str):
        """Write metrics to YAML file"""
        with open(filepath, 'w') as f:
            yaml.dump(metrics, f, default_flow_style=False)

    def run(self):
        """Run the count command"""
        # Validate required arguments
        if not self.args.input_dir:
            raise ValueError("Please provide an input directory")
        if not self.args.output_dir:
            raise ValueError("Please provide an output directory")
        if not self.args.barcode_column:
            raise ValueError("Please provide a barcode column")
        if not self.args.key_column:
            raise ValueError("Please provide a key column")

        # Validate downsampling parameter if provided
        if self.args.downsample is not None:
            if not 0 < self.args.downsample <= 1:
                raise ValueError("Downsample probability must be between 0 and 1")
            if self.args.random_seed is not None:
                random.seed(self.args.random_seed)
                print(f"Using random seed: {self.args.random_seed}", file=sys.stderr)
        
        # Create output directory if it doesn't exist
        os.makedirs(self.args.output_dir, exist_ok=True)
        
        # Read allowed keys if specified
        allowed_keys = None
        if self.args.allowed_list:
            allowed_keys = self._read_allowed_keys(self.args.allowed_list)
            print(f"Loaded {len(allowed_keys)} allowed keys from {self.args.allowed_list}", file=sys.stderr)
        
        # Get list of CSV files in input directory
        input_files = glob.glob(os.path.join(self.args.input_dir, "*.csv"))
        if not input_files:
            raise ValueError(f"No CSV files found in {self.args.input_dir}")
        
        print(f"Found {len(input_files)} CSV files to process", file=sys.stderr)
        
        # Collect metrics for all files
        all_metrics = {}
        
        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(self.args.output_dir, os.path.basename(input_file))
            
            # Process file and raise any errors
            stats = self._process_single_file(
                input_file, output_file, self.args.barcode_column,
                self.args.key_column, self.args.sep, self.args.row_limit,
                allowed_keys, self.args.detailed, self.args.downsample
            )
            
            # Store metrics for this file
            all_metrics[os.path.basename(input_file)] = stats
            
            # Print statistics for this file
            print(f"\nStatistics for {os.path.basename(input_file)}:", file=sys.stderr)
            for category, values in stats.items():
                print(f"\n{category}:", file=sys.stderr)
                if category == 'missing_keys':
                    print(f"  Total missing keys: {values['total_missing']}", file=sys.stderr)
                    if values['total_missing'] > 0 and self.args.detailed:
                        print("  Missing keys:", file=sys.stderr)
                        for key in values['missing_keys_list'][:10]:  # Show first 10 missing keys
                            print(f"    {key}", file=sys.stderr)
                        if len(values['missing_keys_list']) > 10:
                            print(f"    ... and {len(values['missing_keys_list']) - 10} more", file=sys.stderr)
                else:
                    for key, value in values.items():
                        print(f"  {key}: {value:.3f}" if isinstance(value, float) else f"  {key}: {value}", file=sys.stderr)
        
        # Write metrics to YAML file if specified
        if self.args.metrics:
            self._write_metrics(all_metrics, self.args.metrics)
            print(f"\nMetrics written to: {self.args.metrics}", file=sys.stderr)
        
        print(f"\nProcessing complete. Barcode counts written to: {self.args.output_dir}", file=sys.stderr) 