"""Pipeline command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import glob
from pathlib import Path
from typing import List, Tuple
from tqdm import tqdm
import sys

from outerspace.cli.commands.base import BaseCommand
from outerspace.strcmp import get_readpair_files, get_readpairs

class PipelineCommand(BaseCommand):
    """Command for running the complete OUTERSPACE pipeline"""
    
    def __init__(self, args=None):
        from outerspace.cli.main import Cli
        super().__init__(args=args)
        self.Cli = Cli
    
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('pipeline',
            help='Run the complete OUTERSPACE pipeline')
        
        # Required arguments
        parser.add_argument('config_file',
            help='Configuration file with search patterns')
        parser.add_argument('--input-dir',
            help='Directory containing paired FASTQ read files',
            required=True)
        parser.add_argument('--output-dir',
            help='Output directory for all pipeline results',
            required=True)
        
        # Optional arguments
        parser.add_argument('--allowed-list',
            help='Text file containing allowed keys (one per line)')
        parser.add_argument('--mismatches', type=int, default=2,
            help='Number of mismatches allowed for clustering (default: 2)')
        parser.add_argument('--method', choices=['cluster', 'adjacency', 'directional'],
            default='directional',
            help='Clustering method to use (default: directional)')
        parser.add_argument('--barcode-columns', default='UMI_5prime,UMI_3prime',
            help='Column(s) containing barcodes to correct (default: UMI_5prime,UMI_3prime)')
        parser.add_argument('--key-column', default='protospacer',
            help='Column to group by for counting (default: protospacer)')
        parser.add_argument('--sep', default=',',
            help='CSV separator (default: ,)')
        parser.add_argument('--metrics', action='store_true',
            help='Generate metrics files for collapse and count steps')
        parser.add_argument('--config',
            help='TOML configuration file containing command settings')
        return parser

    def _setup_directories(self, output_dir: str) -> Tuple[str, str, str]:
        """Create and return paths for pipeline output directories"""
        extracted_dir = os.path.join(output_dir, 'extracted')
        collapsed_dir = os.path.join(output_dir, 'collapsed')
        counted_dir = os.path.join(output_dir, 'counted')
        
        for directory in [extracted_dir, collapsed_dir, counted_dir]:
            os.makedirs(directory, exist_ok=True)
        
        return extracted_dir, collapsed_dir, counted_dir

    def _run_findseq(self, config_file: str, read1: str, read2: str, output_file: str):
        """Run findseq command for a single pair of reads"""
        findseq_args = [
            'findseq',
            config_file,
            '-1', read1,
            '-2', read2,
            '-o', output_file
        ]
        cli = self.Cli(findseq_args)
        cli.run()

    def _run_collapse(self, input_dir: str, output_dir: str, columns: str, 
                     mismatches: int, method: str, metrics: bool):
        """Run collapse command on all files in input directory"""
        collapse_args = [
            'collapse',
            '--input-dir', input_dir,
            '--output-dir', output_dir,
            '--columns', columns,
            '--mismatches', str(mismatches),
            '--method', method
        ]
        
        if metrics:
            collapse_args.extend(['--metrics', os.path.join(output_dir, 'collapse_metrics.yaml')])
        
        cli = self.Cli(collapse_args)
        cli.run()

    def _run_count(self, input_dir: str, output_dir: str, barcode_column: str,
                  key_column: str, allowed_list: str, metrics: bool):
        """Run count command on all files in input directory"""
        count_args = [
            'count',
            '--input-dir', input_dir,
            '--output-dir', output_dir,
            '--barcode-column', barcode_column,
            '--key-column', key_column
        ]
        
        if allowed_list:
            count_args.extend(['--allowed-list', allowed_list])
        
        if metrics:
            count_args.extend(['--metrics', os.path.join(output_dir, 'count_metrics.yaml')])
        
        cli = self.Cli(count_args)
        cli.run()

    def run(self):
        """Run the pipeline command"""
        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)
        
        # Merge config and args with defaults
        defaults = {
            'mismatches': 2,
            'method': 'directional',
            'barcode_columns': 'UMI_5prime,UMI_3prime',
            'key_column': 'protospacer',
            'sep': ',',
            'metrics': False
        }
        self._merge_config_and_args(defaults)

        # Validate input directory
        if not os.path.exists(self.args.input_dir):
            raise ValueError(f"Input directory not found: {self.args.input_dir}")
        
        # Create output directories
        extracted_dir, collapsed_dir, counted_dir = self._setup_directories(self.args.output_dir)
        
        # Get list of FASTQ files and find read pairs
        fastq_files = glob.glob(os.path.join(self.args.input_dir, "*.fastq.gz"))
        if not fastq_files:
            raise ValueError(f"No FASTQ files found in {self.args.input_dir}")
        
        readpairs = get_readpairs(fastq_files)
        nts = get_readpair_files(readpairs)
        
        print(f"Found {len(nts)} read pairs to process", file=sys.stderr)
        
        # Step 1: Run findseq for each read pair
        for ntd in tqdm(nts, desc="Running findseq"):
            output_file = os.path.join(extracted_dir, f"{ntd.fcsv}.csv")
            self._run_findseq(
                self.args.config_file,
                ntd.fread1,
                ntd.fread2,
                output_file
            )
        
        # Step 2: Run collapse on all extracted files
        print("\nRunning collapse step...", file=sys.stderr)
        self._run_collapse(
            extracted_dir,
            collapsed_dir,
            self.args.barcode_columns,
            self.args.mismatches,
            self.args.method,
            self.args.metrics
        )
        
        # Step 3: Run count on all collapsed files
        print("\nRunning count step...", file=sys.stderr)
        self._run_count(
            collapsed_dir,
            counted_dir,
            f"{self.args.barcode_columns.replace(',', '_')}_corrected",
            self.args.key_column,
            self.args.allowed_list,
            self.args.metrics
        )
        
        print(f"\nPipeline complete. Results written to: {self.args.output_dir}", file=sys.stderr) 