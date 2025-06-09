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
import logging

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.stats import GiniCoefficient
from outerspace.cli.logging_config import setup_logging

class CountCommand(BaseCommand):
    """Command for counting unique barcodes per key value in CSV files"""
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('count',
            help='Count unique barcodes per key value in CSV files')
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument('--input-file',
            help='Input CSV file to process')
        input_group.add_argument('--input-dir',
            help='Input directory containing CSV files to process')
        output_group = parser.add_mutually_exclusive_group(required=True)
        output_group.add_argument('--output-file',
            help='Output CSV file for barcode counts')
        output_group.add_argument('--output-dir',
            help='Output directory for barcode counts')
        parser.add_argument('--barcode-column', required=True,
            help='Column containing barcodes')
        parser.add_argument('--key-column', required=True,
            help='Column to group by')
        parser.add_argument('--sep', default=',',
            help='CSV separator (default: ,)')
        parser.add_argument('--row-limit', type=int,
            help='Process only the first N rows (for testing)',
            default=None)
        parser.add_argument('--allowed-list',
            help='Text file containing allowed keys (one per line)',
            default=None)
        parser.add_argument('--detailed', action='store_true',
            help='Include barcode lists in output')
        parser.add_argument('--downsample', type=float,
            help='Randomly sample reads with probability between 0 and 1')
        parser.add_argument('--random-seed', type=int,
            help='Random seed for downsampling')
        parser.add_argument('--config',
            help='YAML configuration file for command')
        parser.add_argument('--log-file',
            help='Path to log file')
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
        logger = logging.getLogger(__name__)
        
        # Create UMI object for this file
        umi = UMI(mismatches=0)
        key_umi = UMI(mismatches=0)  # For key counts
        
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
                            key_umi.consume(key)
                else:
                    if key and barcode:  # Skip empty values
                        barcodes_by_key[key].add(barcode)
                        umi.consume(barcode)
                        key_umi.consume(key)
        
        # Calculate summary statistics
        total_keys = len(barcodes_by_key)
        total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())
        
        # Calculate Gini coefficients using the stats module
        barcode_result = GiniCoefficient.calculate(umi)
        key_result = GiniCoefficient.calculate(key_umi, allowed_list=list(allowed_keys) if allowed_keys else None)
        barcode_gini = barcode_result
        key_gini = key_result
        
        # Log statistics
        logger.info(f"File statistics for {os.path.basename(input_file)}:")
        logger.info(f"Total rows scanned: {total_rows}")
        logger.info(f"Total keys: {total_keys}")
        logger.info(f"Total barcodes: {total_barcodes}")
        logger.info(f"Average barcodes per key: {total_barcodes / total_keys if total_keys > 0 else 0:.3f}")
        logger.info(f"Barcode Gini coefficient: {barcode_gini:.3f}")
        logger.info(f"Key Gini coefficient: {key_gini:.3f}")
        
        if allowed_keys:
            logger.info(f"Rows with allowed key: {rows_with_allowed_key}")
            missing_keys = allowed_keys - set(barcodes_by_key.keys())
            logger.info(f"Total missing keys: {len(missing_keys)}")
            if len(missing_keys) > 0 and detailed:
                logger.info("Missing keys:")
                for key in sorted(list(missing_keys))[:10]:
                    logger.info(f"  {key}")
                if len(missing_keys) > 10:
                    logger.info(f"  ... and {len(missing_keys) - 10} more")
        
        # Write output
        self._write_counts(barcodes_by_key, output_file, sep, detailed, key_col)
        
        return {
            'total_rows': total_rows,
            'total_keys': total_keys,
            'total_barcodes': total_barcodes,
            'barcode_gini': barcode_gini,
            'key_gini': key_gini
        }

    def _write_counts(self, barcodes_by_key: Dict[str, Set[str]], filepath: str, sep: str, detailed: bool, key_col: str):
        """Write barcode counts per key to CSV file"""
        with open(filepath, 'w', newline='') as f:
            if detailed:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow([key_col, 'unique_barcodes', 'count'])
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, ','.join(sorted(barcodes)), len(barcodes)])
            else:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow([key_col, 'count'])
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, len(barcodes)])

    def run(self):
        """Run the count command"""
        # Set up logging
        logger = setup_logging(log_file=self.args.log_file)
        
        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)
        
        # Merge config and args with defaults
        defaults = {
            'sep': ',',
            'row_limit': None,
            'detailed': False,
            'downsample': None,
            'random_seed': None
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.barcode_column:
            raise ValueError("Please provide a barcode column")
        if not self.args.key_column:
            raise ValueError("Please provide a key column")

        # Validate input/output arguments
        if not self.args.input_file and not self.args.input_dir:
            raise ValueError("Please provide either --input-file or --input-dir")
        if not self.args.output_file and not self.args.output_dir:
            raise ValueError("Please provide either --output-file or --output-dir")

        # Validate downsampling parameter if provided
        if self.args.downsample is not None:
            if not 0 < self.args.downsample <= 1:
                raise ValueError("Downsample probability must be between 0 and 1")
            if self.args.random_seed is not None:
                random.seed(self.args.random_seed)
                logger.info(f"Using random seed: {self.args.random_seed}")
        
        # Read allowed keys if specified
        allowed_keys = None
        if self.args.allowed_list:
            allowed_keys = self._read_allowed_keys(self.args.allowed_list)
            logger.info(f"Loaded {len(allowed_keys)} allowed keys from {self.args.allowed_list}")
        
        # Handle single file case
        if self.args.input_file:
            if not os.path.exists(self.args.input_file):
                raise ValueError(f"Input file not found: {self.args.input_file}")
            
            # Create output directory if needed
            os.makedirs(os.path.dirname(self.args.output_file), exist_ok=True)
            
            try:
                self._process_single_file(
                    self.args.input_file, self.args.output_file,
                    self.args.barcode_column, self.args.key_column,
                    self.args.sep, self.args.row_limit, allowed_keys,
                    self.args.detailed, self.args.downsample
                )
                
            except Exception as e:
                logger.error(f"Error processing {self.args.input_file}: {e}")
                raise e
            
            logger.info(f"Processing complete. Barcode counts written to: {self.args.output_file}")
            return
        
        # Handle directory case
        if not os.path.exists(self.args.input_dir):
            raise ValueError(f"Input directory not found: {self.args.input_dir}")
        
        # Create output directory if it doesn't exist
        os.makedirs(self.args.output_dir, exist_ok=True)
        
        # Get list of CSV files in input directory
        input_files = glob.glob(os.path.join(self.args.input_dir, "*.csv"))
        if not input_files:
            raise ValueError(f"No CSV files found in {self.args.input_dir}")
        
        logger.info(f"Found {len(input_files)} CSV files to process")
        
        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(self.args.output_dir, os.path.basename(input_file))
            
            try:
                self._process_single_file(
                    input_file, output_file, self.args.barcode_column,
                    self.args.key_column, self.args.sep, self.args.row_limit,
                    allowed_keys, self.args.detailed, self.args.downsample
                )
                
            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise e
        
        logger.info(f"Processing complete. Barcode counts written to: {self.args.output_dir}") 