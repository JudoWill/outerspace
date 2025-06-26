"""Barcode correction command for UMI clustering.

This module provides the CollapseCommand class for correcting barcodes in CSV files
using UMI-tools clustering. It supports both single file and batch processing
with various clustering methods and comprehensive metrics reporting.
"""

import csv
import glob
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from tqdm import tqdm

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class CollapseCommand(BaseCommand):
    """Command for correcting barcodes in CSV files using UMI clustering.

    This command performs barcode correction using UMI-tools clustering algorithms.
    It supports processing single files or entire directories, with options for
    different clustering methods and comprehensive metrics reporting.
    """

    def _init_parser(self, subparsers) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparsers
            Subparser group to add command arguments to
        """
        parser = subparsers.add_parser(
            "collapse", help="Correct barcodes in CSV files using UMI-tools clustering"
        )

        # Input options (mutually exclusive)
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--input-file", help="Input CSV file to process")
        input_group.add_argument(
            "--input-dir", help="Input directory containing CSV files to process"
        )

        # Output options (mutually exclusive)
        output_group = parser.add_mutually_exclusive_group(required=True)
        output_group.add_argument(
            "--output-file", help="Output CSV file for corrected barcodes"
        )
        output_group.add_argument(
            "--output-dir", help="Output directory for corrected CSV files"
        )

        # Processing options
        parser.add_argument(
            "--columns",
            help="Column(s) containing barcodes to correct. Can be a single column or comma-separated list",
        )
        parser.add_argument(
            "--mismatches",
            type=int,
            default=2,
            help="Number of mismatches allowed for clustering (default: 2)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--row-limit",
            type=int,
            help="Process only the first N rows (for testing)",
            default=None,
        )
        parser.add_argument(
            "--method",
            choices=["cluster", "adjacency", "directional"],
            default="directional",
            help="Clustering method to use (default: directional)",
        )
        parser.add_argument("--metrics", help="Output YAML file for metrics")
        parser.add_argument(
            "--config", help="TOML configuration file containing command settings"
        )

    def _parse_columns(self, columns_str: str) -> List[str]:
        """Parse comma-separated column string into list of column names.

        Parameters
        ----------
        columns_str : str
            Comma-separated string of column names

        Returns
        -------
        List[str]
            List of column names with whitespace stripped
        """
        return [col.strip() for col in columns_str.split(",")]

    def _process_single_file(
        self,
        input_file: str,
        output_file: str,
        columns: List[str],
        mismatches: int,
        sep: str,
        row_limit: Optional[int],
        method: str,
    ) -> Dict[str, Any]:
        """Process a single CSV file and return metrics.

        This method reads a CSV file, performs barcode correction using UMI clustering,
        and writes the corrected data to an output file.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        output_file : str
            Path to output CSV file
        columns : List[str]
            List of column names containing barcodes
        mismatches : int
            Number of mismatches allowed for clustering
        sep : str
            CSV separator character
        row_limit : Optional[int]
            Maximum number of rows to process (for testing)
        method : str
            Clustering method to use

        Returns
        -------
        Dict[str, Any]
            Dictionary containing processing metrics

        Raises
        ------
        ValueError
            If required columns are not found in the input file
        """
        logger.info(f"Processing file: {input_file}")

        # Read all rows first
        rows = []
        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames

            # Verify columns exist
            missing_cols = [col for col in columns if col not in headers]
            if missing_cols:
                raise ValueError(
                    f"Columns not found in input file: {', '.join(missing_cols)}"
                )

            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break
                rows.append(row)

        if row_limit:
            logger.info(f"Processing first {row_limit} rows of {input_file}")

        # Create UMI object and process barcodes
        umi = UMI(mismatches=mismatches, method=method)

        # Add barcodes to UMI object
        for row in rows:
            # Join multiple columns if specified
            if len(columns) > 1:
                # Check if all columns are present
                if all(row.get(col, "") for col in columns):
                    combined_bc = "".join(str(row[col]) for col in columns)
                    umi.consume(combined_bc)
                    # TODO: maybe there is a way to do "composite" umis
                else:
                    logger.info(f"Skipping row with missing columns: {row}")
            else:
                combined_bc = str(row[columns[0]])
                if combined_bc:  # Skip empty values
                    umi.consume(combined_bc)

        # Create mapping
        logger.info(
            f"Creating clusters from {len(umi._counts)} unique barcodes from {len(rows)} rows "
            f"with {mismatches} mismatches using {method} method"
        )
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
                    corrected_row[key] = corrected.decode("ascii")
                except KeyError:
                    corrected_row[key] = combined_bc
            else:
                corrected_row[key] = ""

            corrected_rows.append(corrected_row)

        # Generate metrics
        metrics = {
            "barcode_counts": {
                "unique_barcodes_before": len(umi._counts),
                "unique_barcodes_after": len(umi.corrected_counts),
                "total_reads": sum(umi._counts.values()),
            },
            "correction_details": {
                "clusters_formed": len(set(umi._mapping.values())),
                "barcodes_corrected": len(umi._mapping) - len(umi.corrected_counts),
            },
        }

        # Write output
        self._write_csv(corrected_rows, output_file, sep)

        logger.info(
            f"Processed {len(rows)} rows, corrected {metrics['correction_details']['barcodes_corrected']} barcodes"
        )
        return metrics

    def _write_csv(self, rows: List[Dict[str, str]], filepath: str, sep: str) -> None:
        """Write rows to CSV file.

        Parameters
        ----------
        rows : List[Dict[str, str]]
            List of dictionaries representing CSV rows
        filepath : str
            Path to output CSV file
        sep : str
            CSV separator character
        """
        if not rows:
            logger.warning("No rows to write")
            return

        logger.debug(f"Writing {len(rows)} rows to {filepath}")

        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter=sep)
            writer.writeheader()
            writer.writerows(rows)

    def _write_metrics(self, metrics: Dict[str, Any], filepath: str) -> None:
        """Write metrics to YAML file.

        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary of metrics to write
        filepath : str
            Path to output YAML file
        """
        import yaml

        logger.info(f"Writing metrics to {filepath}")

        with open(filepath, "w") as f:
            yaml.dump(metrics, f, default_flow_style=False)

    def run(self) -> None:
        """Run the collapse command.

        This method orchestrates the barcode correction process, handling both
        single file and batch processing modes with comprehensive error handling
        and metrics reporting.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        logger.info("Starting barcode correction process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "mismatches": 2,
            "sep": ",",
            "method": "directional",
            "row_limit": None,
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.columns and not self.args.config:
            raise ValueError("Please provide either --columns or --config")

        # Parse columns argument
        columns = self._parse_columns(self.args.columns) if self.args.columns else []

        # Validate input/output arguments
        if not self.args.input_file and not self.args.input_dir:
            raise ValueError("Please provide either --input-file or --input-dir")
        if not self.args.output_file and not self.args.output_dir:
            raise ValueError("Please provide either --output-file or --output-dir")

        # Handle single file case
        if self.args.input_file:
            if not self.args.output_file:
                raise ValueError(
                    "Please provide an output file when using --input-file"
                )

            # Validate input file exists
            if not os.path.exists(self.args.input_file):
                raise ValueError(f"Input file not found: {self.args.input_file}")

            # Create output directory if needed
            output_path = Path(self.args.output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            try:
                metrics = self._process_single_file(
                    self.args.input_file,
                    self.args.output_file,
                    columns,
                    self.args.mismatches,
                    self.args.sep,
                    self.args.row_limit,
                    self.args.method,
                )

                # Print metrics
                logger.info(f"Metrics for {os.path.basename(self.args.input_file)}:")
                for category, values in metrics.items():
                    logger.info(f"{category}:")
                    for key, value in values.items():
                        logger.info(f"  {key}: {value}")

                # Write metrics to YAML file if specified
                if self.args.metrics:
                    self._write_metrics(
                        {os.path.basename(self.args.input_file): metrics},
                        self.args.metrics,
                    )
                    logger.info(f"Metrics written to: {self.args.metrics}")

            except Exception as e:
                logger.error(f"Error processing {self.args.input_file}: {e}")
                raise

            logger.info(
                f"Processing complete. Corrected file written to: {self.args.output_file}"
            )
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

        # Collect metrics for all files
        all_metrics = {}

        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(
                self.args.output_dir, os.path.basename(input_file)
            )

            try:
                metrics = self._process_single_file(
                    input_file,
                    output_file,
                    columns,
                    self.args.mismatches,
                    self.args.sep,
                    self.args.row_limit,
                    self.args.method,
                )

                # Store metrics for this file
                all_metrics[os.path.basename(input_file)] = metrics

                # Print metrics for this file
                logger.info(f"Metrics for {os.path.basename(input_file)}:")
                for category, values in metrics.items():
                    logger.info(f"{category}:")
                    for key, value in values.items():
                        logger.info(f"  {key}: {value}")

            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise

        # Write metrics to YAML file if specified
        if self.args.metrics:
            self._write_metrics(all_metrics, self.args.metrics)
            logger.info(f"Metrics written to: {self.args.metrics}")

        logger.info(
            f"Processing complete. Corrected files written to: {self.args.output_dir}"
        )


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
