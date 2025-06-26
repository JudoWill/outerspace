"""Statistics calculation command for UMI analysis.

This module provides the StatsCommand class for calculating comprehensive statistics
from UMI count data. It supports various diversity metrics, error rates, and
efficiency measures for single-library analysis.
"""

import csv
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from argparse import ArgumentParser

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.stats import (
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    UMIErrorRate,
    UMIRedundancy,
)

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class StatsCommand(BaseCommand):
    """Command for calculating all single-library statistics from counts in a CSV column.

    This command computes comprehensive statistics for UMI count data including
    diversity metrics (Gini coefficient, Shannon diversity, Simpson diversity),
    efficiency measures (recovery rate, efficiency rate, error rate), and
    redundancy analysis.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "stats",
            help="Calculate all single-library statistics from counts in a CSV column",
        )
        parser.add_argument(
            "input_files",
            nargs="+",
            help="Input CSV file(s) to process (supports glob patterns)",
        )
        parser.add_argument("--key-column", help="Column containing keys")
        parser.add_argument(
            "--count-column", help="Column containing pre-counted values"
        )
        parser.add_argument(
            "--scale",
            type=float,
            help="Scale factor for normalized values (e.g., if normalized to mean=1)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--allowed-list", help="Text file containing allowed values (one per line)"
        )
        self._add_common_args(parser)

    def _read_allowed_list(self, filepath: str) -> List[str]:
        """Read allowed values from a text file.

        Parameters
        ----------
        filepath : str
            Path to text file containing allowed values

        Returns
        -------
        List[str]
            List of allowed values with empty lines filtered out

        Raises
        ------
        FileNotFoundError
            If the file doesn't exist
        """
        try:
            with open(filepath, "r") as f:
                allowed_values = [line.strip() for line in f if line.strip()]
            logger.info(f"Loaded {len(allowed_values)} allowed values from {filepath}")
            return allowed_values
        except FileNotFoundError:
            logger.error(f"Allowed list file not found: {filepath}")
            raise

    def _calculate_stats(
        self,
        input_file: str,
        umi_column: str,
        count_column: Optional[str] = None,
        scale: Optional[float] = None,
        sep: str = ",",
        allowed_list: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Calculate all statistics for counts in a CSV column.

        This method computes comprehensive statistics for UMI count data including
        diversity metrics, efficiency measures, and basic count statistics.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        umi_column : str
            Column containing UMIs
        count_column : Optional[str], default=None
            Optional column containing pre-counted values
        scale : Optional[float], default=None
            Optional scale factor for normalized values
        sep : str, default=','
            CSV separator
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed values for filtering

        Returns
        -------
        Dict[str, Any]
            Dictionary containing all calculated statistics

        Notes
        -----
        The method assumes the UMI data is already corrected and focuses on
        statistical analysis rather than clustering.
        """
        logger.info(f"Calculating statistics for {os.path.basename(input_file)}")

        # Create UMI object from CSV
        umi = UMI.from_csv(
            filepath=input_file,
            column=umi_column,
            sep=sep,
            correct=False,  # Assume already corrected
            count_column=count_column,
            scale=scale,
        )

        # Calculate all statistics
        stats = {
            "filename": os.path.basename(input_file),
            "gini_coefficient": GiniCoefficient.calculate(
                umi, allowed_list=allowed_list
            ),
            "shannon_diversity": ShannonDiversity.calculate(
                umi, allowed_list=allowed_list
            ),
            "simpson_diversity": SimpsonDiversity.calculate(
                umi, allowed_list=allowed_list
            ),
            "umi_recovery_rate": UMIRecoveryRate.calculate(
                umi, allowed_list=allowed_list
            ),
            "umi_efficiency_rate": UMIEfficiencyRate.calculate(
                umi, allowed_list=allowed_list
            ),
            "umi_error_rate": UMIErrorRate.calculate(umi),
            "umi_redundancy": UMIRedundancy.calculate(umi, allowed_list=allowed_list),
        }

        # Get some basic statistics for reporting
        counts = umi.corrected_counts
        stats["total_unique"] = len(counts)
        stats["total_counts"] = sum(counts.values())

        # Log all calculated statistics
        logger.info(f"Statistics for {os.path.basename(input_file)}:")
        for stat_name, value in stats.items():
            if stat_name not in ["filename", "total_unique", "total_counts"]:
                if value is not None:
                    logger.info(f"  {stat_name}: {value:.4f}")
                else:
                    logger.info(f"  {stat_name}: None")

        if allowed_list:
            missing = set(allowed_list) - {k.decode("ascii") for k in counts.keys()}
            logger.info(f"  Missing values: {len(missing)}/{len(allowed_list)}")
            if missing and len(missing) <= 10:
                logger.info("  Missing values list:")
                for value in sorted(missing):
                    logger.info(f"    {value}")
            elif missing:
                logger.info(
                    f"  Missing values list: {len(missing)} values (showing first 10)"
                )
                for value in sorted(list(missing)[:10]):
                    logger.info(f"    {value}")

        return stats

    def run(self) -> None:
        """Run the stats command.

        This method orchestrates the statistics calculation process, handling
        configuration loading, file processing, and result output with comprehensive
        error handling and logging.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid, or if no files are
            successfully processed
        """
        logger.info("Starting statistics calculation process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {"sep": ",", "scale": None}
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.key_column and not self.args.config:
            raise ValueError("Please provide either --key-column or --config")

        # Validate scale if provided
        if self.args.scale is not None and self.args.scale <= 0:
            raise ValueError("Scale value must be positive")

        try:
            # Read allowed list if provided
            allowed_list = None
            if self.args.allowed_list:
                allowed_list = self._read_allowed_list(self.args.allowed_list)

            # Process each input file
            all_stats = []
            processed_files = 0

            for input_file in self.args.input_files:
                if not os.path.exists(input_file):
                    logger.warning(f"Input file not found: {input_file}")
                    continue

                try:
                    stats = self._calculate_stats(
                        input_file,
                        self.args.key_column,
                        count_column=self.args.count_column,
                        scale=self.args.scale,
                        sep=self.args.sep,
                        allowed_list=allowed_list,
                    )
                    all_stats.append(stats)
                    processed_files += 1

                except Exception as e:
                    logger.error(f"Error processing {input_file}: {e}")
                    continue

            if not all_stats:
                raise ValueError("No files were successfully processed")

            logger.info(f"Successfully processed {processed_files} files")

            # Write all results to stdout as CSV
            writer = csv.DictWriter(sys.stdout, fieldnames=all_stats[0].keys())
            writer.writeheader()
            writer.writerows(all_stats)

            logger.info("Statistics calculation completed successfully")

        except Exception as e:
            logger.error(f"Error calculating statistics: {e}")
            raise ValueError(f"Error calculating statistics: {e}")


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
