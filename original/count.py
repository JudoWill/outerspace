#!/usr/bin/env python3
"""Script to count unique barcodes per key value in CSV files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import argparse
import csv
from collections import Counter, defaultdict
from typing import List, Dict, Any, Set
import sys
from tqdm import tqdm
import os
import glob

def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Count unique barcodes per key value in CSV files")
    parser.add_argument("input_dir", help="Input directory containing CSV files")
    parser.add_argument("output_dir", help="Output directory for barcode counts")
    parser.add_argument("--barcode-column", required=True, help="Column containing barcodes")
    parser.add_argument("--key-column", required=True, help="Column to group by")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--row-limit", type=int, help="Process only the first N rows (for testing)")
    parser.add_argument("--allowed-list", help="Text file containing allowed keys (one per line)")
    parser.add_argument("--detailed", action="store_true", help="Include barcode lists in output (default: False)")
    return parser.parse_args()

def read_allowed_keys(filepath: str) -> Set[str]:
    """Read allowed keys from a text file"""
    allowed_keys = set()
    with open(filepath, 'r') as f:
        for line in f:
            key = line.strip()
            if key:  # Skip empty lines
                allowed_keys.add(key)
    return allowed_keys

def get_barcodes_per_key(filepath: str, barcode_col: str, key_col: str, sep: str, row_limit: int = None, allowed_keys: Set[str] = None) -> Dict[str, Set[str]]:
    """Get unique barcodes for each key value"""
    barcodes_by_key = defaultdict(set)
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter=sep)
        headers = reader.fieldnames
        
        # Verify columns exist
        missing_cols = [col for col in [barcode_col, key_col] if col not in headers]
        if missing_cols:
            raise ValueError(f"Columns not found in input file: {', '.join(missing_cols)}")
        
        for i, row in enumerate(tqdm(reader, desc="Reading rows")):
            if row_limit and i >= row_limit:
                break
                
            key = str(row[key_col])
            barcode = str(row[barcode_col])
            
            # Skip if key is not in allowed list
            if allowed_keys and key not in allowed_keys:
                continue
                
            if key and barcode:  # Skip empty values
                barcodes_by_key[key].add(barcode)
    
    return barcodes_by_key

def write_counts(barcodes_by_key: Dict[str, Set[str]], filepath: str, sep: str, detailed: bool = False):
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

def process_single_file(input_file: str, output_file: str, barcode_col: str, key_col: str,
                       sep: str, row_limit: int, allowed_keys: Set[str], detailed: bool) -> Dict[str, Any]:
    """Process a single CSV file and return summary statistics"""
    # Get barcode counts per key
    if row_limit:
        print(f"Processing first {row_limit} rows of {input_file}", file=sys.stderr)
    
    barcodes_by_key = get_barcodes_per_key(
        input_file, 
        barcode_col, 
        key_col, 
        sep, 
        row_limit,
        allowed_keys
    )
    
    # Calculate summary statistics
    total_keys = len(barcodes_by_key)
    total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())
    
    # Write output
    write_counts(barcodes_by_key, output_file, sep, detailed)
    
    return {
        'file_stats': {
            'total_keys': total_keys,
            'total_barcodes': total_barcodes,
            'average_barcodes_per_key': total_barcodes / total_keys if total_keys > 0 else 0
        }
    }

def main():
    args = parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Read allowed keys if specified
        allowed_keys = None
        if args.allowed_list:
            allowed_keys = read_allowed_keys(args.allowed_list)
            print(f"Loaded {len(allowed_keys)} allowed keys from {args.allowed_list}", file=sys.stderr)
        
        # Get list of CSV files in input directory
        input_files = glob.glob(os.path.join(args.input_dir, "*.csv"))
        if not input_files:
            print(f"No CSV files found in {args.input_dir}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(input_files)} CSV files to process", file=sys.stderr)
        
        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(args.output_dir, os.path.basename(input_file))
            
            try:
                stats = process_single_file(
                    input_file, output_file, args.barcode_column,
                    args.key_column, args.sep, args.row_limit,
                    allowed_keys, args.detailed
                )
                
                # Print statistics for this file
                print(f"\nStatistics for {os.path.basename(input_file)}:", file=sys.stderr)
                for category, values in stats.items():
                    print(f"\n{category}:", file=sys.stderr)
                    for key, value in values.items():
                        print(f"  {key}: {value}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error processing {input_file}: {e}", file=sys.stderr)
                continue
        
        print(f"\nProcessing complete. Barcode counts written to: {args.output_dir}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 