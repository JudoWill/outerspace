"""Statistics calculation command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import csv
import os
import sys
import logging
from typing import List, Optional, Dict, Any
from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.stats import (
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    UMIErrorRate,
    UMIRedundancy
)

class StatsCommand(BaseCommand):
    """Command for calculating all single-library statistics from counts in a CSV column"""
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('stats',
            help='Calculate all single-library statistics from counts in a CSV column')
        parser.add_argument('input_files', nargs='+',
            help='Input CSV file(s) to process (supports glob patterns)')
        parser.add_argument('--umi-column', required=True,
            help='Column containing UMIs')
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

    def _calculate_stats(self, input_file: str, umi_column: str, count_column: Optional[str] = None, 
                        scale: Optional[float] = None, sep: str = ",", 
                        allowed_list: Optional[List[str]] = None) -> Dict[str, Any]:
        """Calculate all statistics for counts in a CSV column
        
        Args:
            input_file: Path to input CSV file
            umi_column: Column containing UMIs
            count_column: Optional column containing pre-counted values
            scale: Optional scale factor for normalized values
            sep: CSV separator
            allowed_list: Optional list of allowed values
            
        Returns:
            Dictionary containing all calculated statistics
        """
        logger = logging.getLogger(__name__)
        
        # Create UMI object from CSV
        umi = UMI.from_csv(
            filepath=input_file,
            column=umi_column,
            sep=sep,
            correct=False,  # Assume already corrected
            count_column=count_column,
            scale=scale
        )
        
        # Calculate all statistics
        stats = {
            'filename': os.path.basename(input_file),
            'gini_coefficient': GiniCoefficient.calculate(umi, allowed_list=allowed_list),
            'shannon_diversity': ShannonDiversity.calculate(umi, allowed_list=allowed_list),
            'simpson_diversity': SimpsonDiversity.calculate(umi, allowed_list=allowed_list),
            'umi_recovery_rate': UMIRecoveryRate.calculate(umi, allowed_list=allowed_list),
            'umi_efficiency_rate': UMIEfficiencyRate.calculate(umi, allowed_list=allowed_list),
            'umi_error_rate': UMIErrorRate.calculate(umi),
            'umi_redundancy': UMIRedundancy.calculate(umi, allowed_list=allowed_list)
        }
        
        # Get some basic statistics for reporting
        counts = umi.corrected_counts
        stats['total_unique'] = len(counts)
        stats['total_counts'] = sum(counts.values())
                
        # Log all calculated statistics
        for stat_name, value in stats.items():
            if stat_name not in ['filename', 'total_unique', 'total_counts']:
                logger.info(f"  {stat_name}: {value:.4f}" if value is not None else f"  {stat_name}: None")
        
        if allowed_list:
            missing = set(allowed_list) - {k.decode('ascii') for k in counts.keys()}
            logger.info(f"  Missing values: {len(missing)}/{len(allowed_list)}")
            if missing and len(missing) <= 10:
                logger.info("  Missing values list:")
                for value in sorted(missing):
                    logger.info(f"    {value}")
        
        return stats

    def run(self):
        """Run the stats command"""
        logger = logging.getLogger(__name__)
        
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
        if not self.args.umi_column:
            raise ValueError("Please provide a UMI column")

        # Validate scale if provided
        if self.args.scale is not None and self.args.scale <= 0:
            raise ValueError("Scale value must be positive")

        try:
            # Read allowed list if provided
            allowed_list = None
            if self.args.allowed_list:
                allowed_list = self._read_allowed_list(self.args.allowed_list)
                logger.info(f"Loaded {len(allowed_list)} allowed values from {self.args.allowed_list}")
            
            # Process each input file
            all_stats = []
            for input_file in self.args.input_files:
                if not os.path.exists(input_file):
                    logger.warning(f"Input file not found: {input_file}")
                    continue
                
                try:
                    stats = self._calculate_stats(
                        input_file,
                        self.args.umi_column,
                        count_column=self.args.count_column,
                        scale=self.args.scale,
                        sep=self.args.sep,
                        allowed_list=allowed_list
                    )
                    all_stats.append(stats)
                except Exception as e:
                    logger.error(f"Error processing {input_file}: {e}")
                    continue
            
            if not all_stats:
                raise ValueError("No files were successfully processed")
            
            # Write all results to stdout as CSV
            writer = csv.DictWriter(sys.stdout, fieldnames=all_stats[0].keys())
            writer.writeheader()
            writer.writerows(all_stats)
                
        except Exception as e:
            logger.error(f"Error calculating statistics: {e}")
            raise ValueError(f"Error calculating statistics: {e}") 