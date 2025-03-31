#!/usr/bin/env python3
"""Script to visualize barcode counts from CSV files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Dict, Any
import sys
import os
import glob
from tqdm import tqdm

def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Visualize barcode counts from CSV files")
    parser.add_argument("input_dir", help="Input directory containing CSV files with barcode counts")
    parser.add_argument("output_dir", help="Output directory for visualization plots")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--bins", type=int, default=50, help="Number of histogram bins (default: 50)")
    parser.add_argument("--title-prefix", help="Prefix for plot titles (default: filename)")
    parser.add_argument("--xlabel", default="Number of Unique Barcodes", help="X-axis label")
    parser.add_argument("--ylabel", default="Count", help="Y-axis label")
    parser.add_argument("--log-scale", action="store_true", help="Use log scale for y-axis")
    parser.add_argument("--format", default="png", help="Output image format (default: png)")
    return parser.parse_args()

def read_counts(filepath: str, sep: str) -> List[int]:
    """Read barcode counts from CSV file"""
    counts = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            counts.append(int(row['barcode_count']))
    return counts

def create_histogram(counts: List[int], args: argparse.Namespace, title: str) -> Tuple[plt.Figure, plt.Axes]:
    """Create histogram of barcode counts"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create histogram
    ax.hist(counts, bins=args.bins, edgecolor='black')
    
    # Customize plot
    ax.set_title(title)
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    stats_text = f"Mean: {np.mean(counts):.1f}\n"
    stats_text += f"Median: {np.median(counts):.1f}\n"
    stats_text += f"Max: {np.max(counts)}\n"
    stats_text += f"Total samples: {len(counts)}"
    
    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Set log scale if requested
    if args.log_scale:
        ax.set_yscale('log')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig, ax

def process_single_file(input_file: str, output_file: str, args: argparse.Namespace) -> Dict[str, Any]:
    """Process a single CSV file and create visualization"""
    # Read counts
    counts = read_counts(input_file, args.sep)
    
    if not counts:
        raise ValueError("No counts found in input file")
    
    # Create title
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    title = f"{args.title_prefix} - {base_name}" if args.title_prefix else base_name
    
    # Create histogram
    fig, ax = create_histogram(counts, args, title)
    
    # Save plot
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return {
        'file_stats': {
            'total_samples': len(counts),
            'mean_count': np.mean(counts),
            'median_count': np.median(counts),
            'max_count': np.max(counts),
            'min_count': np.min(counts)
        }
    }

def main():
    args = parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Get list of CSV files in input directory
        input_files = glob.glob(os.path.join(args.input_dir, "*.csv"))
        if not input_files:
            print(f"No CSV files found in {args.input_dir}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(input_files)} CSV files to process", file=sys.stderr)
        
        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            output_file = os.path.join(args.output_dir, f"{base_name}.{args.format}")
            
            try:
                stats = process_single_file(input_file, output_file, args)
                
                # Print statistics for this file
                print(f"\nStatistics for {os.path.basename(input_file)}:", file=sys.stderr)
                for category, values in stats.items():
                    print(f"\n{category}:", file=sys.stderr)
                    for key, value in values.items():
                        print(f"  {key}: {value:.1f}" if isinstance(value, float) else f"  {key}: {value}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error processing {input_file}: {e}", file=sys.stderr)
                continue
        
        print(f"\nProcessing complete. Plots written to: {args.output_dir}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 