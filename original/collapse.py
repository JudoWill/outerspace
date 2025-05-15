#!/usr/bin/env python3
"""Script to correct barcodes in CSV files using UMI-tools clustering"""

__author__ = "Will Dampier"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import argparse
import csv
from collections import Counter
from typing import List, Dict, Any, Tuple
import sys
from tqdm import tqdm
import os
import glob
from outerspace.umi import UMI

def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Correct barcodes in CSV files using UMI-tools clustering")
    parser.add_argument("input_dir", help="Input directory containing CSV files")
    parser.add_argument("output_dir", help="Output directory for corrected CSV files")
    parser.add_argument("--columns", required=True, help="Column(s) containing barcodes to correct. Can be a single column or comma-separated list")
    parser.add_argument("--mismatches", type=int, default=2, help="Number of mismatches allowed for clustering")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--row-limit", type=int, help="Process only the first N rows (for testing)")
    parser.add_argument("--method", choices=["cluster", "adjacency", "directional"], 
                       default="adjacency", help="Clustering method to use (default: directional)")
    return parser.parse_args()

def parse_columns(columns_str: str) -> List[str]:
    """Parse comma-separated column string into list of column names"""
    return [col.strip() for col in columns_str.split(",")]

def process_single_file(input_file: str, output_file: str, columns: List[str], 
                       mismatches: int, sep: str, row_limit: int, method: str) -> Dict[str, Any]:
    """Process a single CSV file and return metrics"""
    # Read all rows first
    rows = []
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=sep)
        headers = reader.fieldnames
        
        # Verify columns exist
        missing_cols = [col for col in columns if col not in headers]
        if missing_cols:
            raise ValueError(f"Columns not found in input file: {', '.join(missing_cols)}")
        
        for i, row in enumerate(tqdm(reader, desc="Reading rows")):
            if row_limit and i >= row_limit:
                break
            rows.append(row)
    
    if row_limit:
        print(f"Processing first {row_limit} rows of {input_file}", file=sys.stderr)
    
    # Create UMI object and process barcodes
    umi = UMI(mismatches=mismatches, method=method)
    
    # Add barcodes to UMI object
    for row in rows:
        # Join multiple columns if specified
        if len(columns) > 1:
            combined_bc = "".join(str(row[col]) for col in columns)
        else:
            combined_bc = str(row[columns[0]])
        if combined_bc:  # Skip empty values
            umi.consume(combined_bc)
    
    # Create mapping
    print(f"Creating clusters from {len(umi._counts)} unique barcodes from {len(rows)} rows with {mismatches} mismatches using {method} method", file=sys.stderr)
    umi.create_mapping()
    
    # Correct barcodes in rows
    corrected_rows = []
    key = "_".join(columns) + "_corrected"
    for row in tqdm(rows, desc="Correcting barcodes"):
        corrected_row = row.copy()
        
        # Join multiple columns if specified
        if len(columns) > 1:
            combined_bc = "".join(str(row[col]) for col in columns)
        else:
            combined_bc = str(row[columns[0]])
        
        # Create corrected column
        if combined_bc:
            try:
                corrected = umi[combined_bc]
                corrected_row[key] = corrected.decode('ascii')
            except KeyError:
                corrected_row[key] = combined_bc
        else:
            corrected_row[key] = ""
        
        corrected_rows.append(corrected_row)
    
    # Generate metrics
    metrics = {
        'barcode_counts': {
            'unique_barcodes_before': len(umi._counts),
            'unique_barcodes_after': len(umi.corrected_counts),
            'total_reads': sum(umi._counts.values())
        },
        'correction_details': {
            'clusters_formed': len(set(umi._mapping.values())),
            'barcodes_corrected': len(umi._mapping) - len(umi.corrected_counts)
        }
    }
    
    # Write output
    write_csv(corrected_rows, output_file, sep)
    
    return metrics

def write_csv(rows: List[Dict[str, str]], filepath: str, sep: str):
    """Write rows to CSV file"""
    if not rows:
        return
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter=sep)
        writer.writeheader()
        writer.writerows(rows)

def main():
    args = parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Parse columns argument
        columns = parse_columns(args.columns)
        
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
                metrics = process_single_file(
                    input_file, output_file, columns, args.mismatches,
                    args.sep, args.row_limit, args.method
                )
                
                # Print metrics for this file
                print(f"\nMetrics for {os.path.basename(input_file)}:", file=sys.stderr)
                for category, values in metrics.items():
                    print(f"\n{category}:", file=sys.stderr)
                    for key, value in values.items():
                        print(f"  {key}: {value}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error processing {input_file}: {e}", file=sys.stderr)
                continue
        
        print(f"\nProcessing complete. Corrected files written to: {args.output_dir}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 