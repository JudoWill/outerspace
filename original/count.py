#!/usr/bin/env python3
"""Script to count unique barcodes per key value in a CSV file"""

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

def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Count unique barcodes per key value in a CSV file")
    parser.add_argument("input", help="Input CSV file containing barcodes")
    parser.add_argument("output", help="Output CSV file for barcode counts")
    parser.add_argument("--barcode-column", required=True, help="Column containing barcodes")
    parser.add_argument("--key-column", required=True, help="Column to group by")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--row-limit", type=int, help="Process only the first N rows (for testing)")
    return parser.parse_args()

def get_barcodes_per_key(filepath: str, barcode_col: str, key_col: str, sep: str, row_limit: int = None) -> Dict[str, Set[str]]:
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
            
            if key and barcode:  # Skip empty values
                barcodes_by_key[key].add(barcode)
    
    return barcodes_by_key

def write_counts(barcodes_by_key: Dict[str, Set[str]], filepath: str, sep: str):
    """Write barcode counts per key to CSV file"""
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow(['key', 'unique_barcodes', 'barcode_count'])
        for key, barcodes in sorted(barcodes_by_key.items()):
            writer.writerow([key, ','.join(sorted(barcodes)), len(barcodes)])

def main():
    args = parse_args()
    
    try:
        # Get barcode counts per key
        if args.row_limit:
            print(f"Processing first {args.row_limit} rows", file=sys.stderr)
        
        barcodes_by_key = get_barcodes_per_key(
            args.input, 
            args.barcode_column, 
            args.key_column, 
            args.sep, 
            args.row_limit
        )
        
        # Print summary
        total_keys = len(barcodes_by_key)
        total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())
        print(f"\nFound {total_keys} unique keys", file=sys.stderr)
        print(f"Total unique barcodes across all keys: {total_barcodes}", file=sys.stderr)
        
        # Write output
        write_counts(barcodes_by_key, args.output, args.sep)
        print(f"\nBarcode counts written to: {args.output}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 