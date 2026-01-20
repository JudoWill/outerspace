"""UMI alignment command using spoa.

This module provides the AlignCommand class for aligning UMI sequences from CSV files
using spoa. It supports both single file and batch processing with configurable
alignment parameters.
"""

import csv
import glob
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from argparse import ArgumentParser

from tqdm import tqdm

from outerspace.cli.commands.base import BaseCommand
from outerspace.align import align_umi_sequences

# Set up logging
logger = logging.getLogger(__name__)

# Increase CSV field size limit to handle large fields
csv.field_size_limit(sys.maxsize)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class AlignCommand(BaseCommand):
    """Command for aligning UMI sequences from CSV files using spoa.

    This command reads UMI sequences from a CSV column, aligns them using spoa,
    and outputs the aligned sequences. It supports processing single files or
    entire directories with configurable alignment parameters.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "align", help="Align UMI sequences from CSV files using spoa"
        )

        # Input options (mutually exclusive)
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--input-file", help="Input CSV file to process")
        input_group.add_argument(
            "--input-dir", help="Input directory containing CSV files to process"
        )

        # Output options (optional - defaults to stdout)
        parser.add_argument(
            "--output-file",
            help="Output file for aligned sequences (default: stdout)",
            default=None,
        )
        parser.add_argument(
            "--output-dir",
            help="Output directory for aligned sequences (for batch processing)",
            default=None,
        )

        # Processing options
        parser.add_argument(
            "--column",
            help="Column name containing UMI sequences to align",
            required=True,
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--row-limit",
            type=int,
            help="Process only the first N rows (for testing)",
            default=None,
        )
        parser.add_argument(
            "--unique",
            action="store_true",
            help="Remove duplicate sequences",
        )


        # Alignment parameters
        parser.add_argument(
            "--match",
            type=int,
            default=5,
            help="Score for matching bases (default: 5)",
        )
        parser.add_argument(
            "--mismatch",
            type=int,
            default=-4,
            help="Penalty for mismatching bases (default: -4)",
        )
        parser.add_argument(
            "--gap",
            type=int,
            default=-8,
            help="Penalty for gaps/indels (default: -8)",
        )
        parser.add_argument(
            "--algorithm",
            type=int,
            choices=[0, 1, 2],
            default=1,
            help=(
                "Alignment algorithm: 0=local (Smith-Waterman), "
                "1=global (Needleman-Wunsch), 2=semi-global (default: 1)"
            ),
        )

        self._add_common_args(parser)

    def _process_single_file(
        self,
        input_file: str,
        output_file: Optional[str],
        column: str,
        sep: str,
        row_limit: Optional[int],
        match: int,
        mismatch: int,
        gap: int,
        algorithm: int,
        unique: bool,
    ) -> Dict[str, Any]:
        """Process a single CSV file and align UMIs.

        This method reads a CSV file, extracts UMI sequences from the specified
        column, aligns them using spoa, and writes the aligned sequences to
        the output file or stdout.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        output_file : Optional[str]
            Path to output file (None for stdout)
        column : str
            Column name containing UMIs
        sep : str
            CSV separator character
        row_limit : Optional[int]
            Maximum number of rows to process (for testing)
        match : int
            Match score for alignment
        mismatch : int
            Mismatch penalty for alignment
        gap : int
            Gap penalty for alignment
        algorithm : int
            Alignment algorithm (0=local, 1=global, 2=semi-global)
        unique : bool
            Remove duplicate sequences
        Returns
        -------
        Dict[str, Any]
            Dictionary containing processing metrics

        Raises
        ------
        ValueError
            If required column is not found in the input file
        """
        logger.info(f"Processing file: {input_file}")

        # Read UMI sequences from CSV
        sequences = []
        unique_sequences = set()
        total_rows = 0

        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames

            # Verify column exists
            if column not in headers:
                raise ValueError(
                    f"Column '{column}' not found in input file. "
                    f"Available columns: {', '.join(headers)}"
                )
            
            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break

                total_rows += 1
                umi_seq = str(row.get(column, "")).strip()
                if umi_seq:  # Only add non-empty sequences
                    if unique:
                        if umi_seq not in unique_sequences:
                            sequences.append(umi_seq)
                            unique_sequences.add(umi_seq)
                    else:
                        sequences.append(umi_seq)

        logger.info(f"Read {len(sequences)} UMI sequences from {total_rows} rows")

        if not sequences:
            logger.warning("No UMI sequences found in column")
            aligned_sequences = []
        else:
            # Perform alignment
            logger.info(
                f"Aligning {len(sequences)} sequences with "
                f"match={match}, mismatch={mismatch}, gap={gap}, algorithm={algorithm}"
            )
            aligned_sequences = align_umi_sequences(
                sequences=sequences,
                match=match,
                mismatch=mismatch,
                gap=gap,
                algorithm=algorithm,
            )

        # Write output
        self._write_aligned_sequences(aligned_sequences, output_file)

        metrics = {
            "total_rows": total_rows,
            "sequences_read": len(sequences),
            "aligned_sequences": len(aligned_sequences),
        }

        if aligned_sequences:
            metrics["aligned_length"] = len(aligned_sequences[0])

        logger.info(f"Processed {len(sequences)} sequences, output {len(aligned_sequences)} aligned sequences")
        return metrics

    def _write_aligned_sequences(
        self, aligned_sequences: List[str], output_file: Optional[str]
    ) -> None:
        """Write aligned sequences to file or stdout.

        Parameters
        ----------
        aligned_sequences : List[str]
            List of aligned sequences to write
        output_file : Optional[str]
            Path to output file (None for stdout)
        """
        if output_file:
            logger.info(f"Writing aligned sequences to {output_file}")
            # Create output directory if needed
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, "w") as f:
                for seq in aligned_sequences:
                    f.write(seq + "\n")
        else:
            # Write to stdout
            logger.debug("Writing aligned sequences to stdout")
            for seq in aligned_sequences:
                print(seq)

    def run(self) -> None:
        """Run the align command.

        This method orchestrates the alignment process, handling both
        single file and batch processing modes with comprehensive error handling.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        logger.info("Starting UMI alignment process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "sep": ",",
            "row_limit": None,
            "match": 5,
            "mismatch": -4,
            "gap": -8,
            "algorithm": 1,
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.column:
            raise ValueError("Please provide --column")

        # Validate input/output arguments
        if not self.args.input_file and not self.args.input_dir:
            raise ValueError("Please provide either --input-file or --input-dir")

        # Handle single file case
        if self.args.input_file:
            if not os.path.exists(self.args.input_file):
                raise ValueError(f"Input file not found: {self.args.input_file}")

            # For single file, output_file is optional (defaults to stdout)
            # but output_dir doesn't make sense
            if self.args.output_dir:
                raise ValueError(
                    "Cannot use --output-dir with --input-file. "
                    "Use --output-file or omit for stdout."
                )

            try:
                metrics = self._process_single_file(
                    input_file=self.args.input_file,
                    output_file=self.args.output_file,
                    column=self.args.column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    match=self.args.match,
                    mismatch=self.args.mismatch,
                    gap=self.args.gap,
                    algorithm=self.args.algorithm,
                    unique=self.args.unique,
                )

                logger.info(f"Alignment complete. Processed {metrics['sequences_read']} sequences")

            except Exception as e:
                logger.error(f"Error processing {self.args.input_file}: {e}")
                raise

            return

        # Handle directory case
        if not self.args.output_dir:
            raise ValueError(
                "Please provide --output-dir when using --input-dir for batch processing"
            )

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
            output_file = os.path.join(
                self.args.output_dir, os.path.basename(input_file)
            )

            try:
                metrics = self._process_single_file(
                    input_file=input_file,
                    output_file=output_file,
                    column=self.args.column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    match=self.args.match,
                    mismatch=self.args.mismatch,
                    gap=self.args.gap,
                    algorithm=self.args.algorithm,
                )

                logger.info(
                    f"Processed {os.path.basename(input_file)}: "
                    f"{metrics['sequences_read']} sequences aligned"
                )

            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise

        logger.info(
            f"Processing complete. Aligned sequences written to: {self.args.output_dir}"
        )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
