#!/usr/bin/env python3
"""Script to calculate Gini coefficient from counts in a CSV column"""

__author__ = "Will Dampier"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import argparse
import csv
import sys
from grna_extraction.umi import UMI
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Calculate Gini coefficient from counts in a CSV column")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("--column", required=True, help="Column to calculate Gini coefficient from")
    parser.add_argument("--count-column", help="Column containing pre-counted values")
    parser.add_argument("--scale", type=float, help="Scale factor for normalized values (e.g., if normalized to mean=1)")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--allowed-list", help="Text file containing allowed values (one per line)")
    return parser.parse_args()


def read_allowed_list(filepath: str) -> list[str]:
    """Read allowed values from a text file"""
    with open(filepath, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def calculate_gini(input_file: str, column: str, count_column: str = None, scale: float = None,
                  sep: str = ",", allowed_list: list[str] = None) -> float:
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
    # Create UMI object (with no correction)
    umi = UMI(mismatches=0)
    
    # Read CSV and count values
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=sep)
        
        # Verify columns exist
        if column not in reader.fieldnames:
            raise ValueError(f"Column '{column}' not found in CSV file")
        if count_column and count_column not in reader.fieldnames:
            raise ValueError(f"Count column '{count_column}' not found in CSV file")
        
        for row in reader:
            value = row[column]
            if not value:  # Skip empty values
                continue
                
            if count_column:
                try:
                    count = float(row[count_column])
                    if scale:
                        count = int(round(count * scale))
                    else:
                        count = int(round(count))
                except (ValueError, TypeError):
                    print(f"Warning: Invalid count value '{row[count_column]}' for {value}", file=sys.stderr)
                    continue
                
                # Add the value count times
                for _ in range(count):
                    umi.consume(value)
            else:
                umi.consume(value)
    
    # Calculate Gini coefficient
    gini = umi.gini_coefficient(allowed_list=allowed_list)
    
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


def main():
    """Main entry point"""
    args = parse_args()
    
    try:
        # Read allowed list if provided
        allowed_list = None
        if args.allowed_list:
            allowed_list = read_allowed_list(args.allowed_list)
            print(f"Loaded {len(allowed_list)} allowed values from {args.allowed_list}", file=sys.stderr)
        
        # Calculate Gini coefficient
        gini = calculate_gini(
            args.input_file,
            args.column,
            count_column=args.count_column,
            scale=args.scale,
            sep=args.sep,
            allowed_list=allowed_list
        )
        
        # Print result to stdout (for scripting)
        if gini is not None:
            print(f"{gini:.6f}")
        else:
            print("NA")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 