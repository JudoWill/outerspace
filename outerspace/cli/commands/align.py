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
from collections import defaultdict
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
    """Command for aligning sequences from CSV files using spoa.

    This command counts unique barcodes per key, filters keys based on barcode count,
    and aligns sequences using spoa. It supports processing single files or entire
    directories with configurable alignment parameters and filtering options.
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
            "--key-column",
            help="Column containing sequences to align",
            default=None,
        )
        parser.add_argument(
            "--barcode-column",
            help="Column containing unique markers for counting",
            default=None,
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--row-limit",
            type=int,
            help="Process only the first N rows (for testing)",
            default=None,
        )
        parser.add_argument(
            "--min-count",
            type=int,
            default=0,
            help="Minimum unique barcode count threshold for a key to be included (default: 0)",
        )
        parser.add_argument(
            "--top-n",
            type=int,
            default=None,
            help="Keep only top N keys by unique barcode count (default: all)",
        )
        parser.add_argument(
            "--min-frequency",
            type=float,
            default=0.0,
            help="Minimum frequency percentage threshold based on unique barcode count (default: 0.0)",
        )
        parser.add_argument(
            "--align-by-barcode",
            action="store_true",
            help="If set, align keys separately grouped by identical barcodes (default: align all together)",
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
        key_col: str,
        barcode_col: str,
        sep: str,
        row_limit: Optional[int],
        match: int,
        mismatch: int,
        gap: int,
        algorithm: int,
        min_count: int,
        top_n: Optional[int],
        min_frequency: float,
        align_by_barcode: bool,
    ) -> Dict[str, Any]:
        """Process a single CSV file and align UMIs.

        This method reads a CSV file, counts unique barcodes per key, filters keys
        based on barcode count, and aligns the sequences using spoa.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        output_file : Optional[str]
            Path to output file (None for stdout)
        key_col : str
            Column name containing sequences to align
        barcode_col : str
            Column name containing unique markers for counting
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
        min_count : int
            Minimum unique barcode count threshold
        top_n : Optional[int]
            Keep only top N keys by unique barcode count
        min_frequency : float
            Minimum frequency percentage threshold
        align_by_barcode : bool
            If True, align keys separately grouped by barcode
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

        # Count unique barcodes per key
        keys_to_barcodes = defaultdict(set)  # key -> set of unique barcodes
        keys_to_sequences = defaultdict(list)  # key -> list of sequences (for alignment)
        total_rows = 0

        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames

            # Verify columns exist
            missing_cols = [col for col in [key_col, barcode_col] if col not in headers]
            if missing_cols:
                raise ValueError(
                    f"Columns not found in input file: {', '.join(missing_cols)}. "
                    f"Available columns: {', '.join(headers)}"
                )
            
            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break

                total_rows += 1
                key = str(row.get(key_col, "")).strip()
                barcode = str(row.get(barcode_col, "")).strip()
                
                if not key or not barcode:
                    continue
                
                # The key itself is the sequence to align
                keys_to_barcodes[key].add(barcode)
                keys_to_sequences[key].append(key)

        logger.info(f"Read {total_rows} rows, found {len(keys_to_barcodes)} unique keys")

        # Calculate unique barcode count per key
        key_counts = {key: len(barcodes) for key, barcodes in keys_to_barcodes.items()}
        
        # Apply filtering
        filtered_keys = self._filter_keys(
            key_counts, min_count, top_n, min_frequency
        )
        
        logger.info(
            f"Filtered {len(filtered_keys)} keys from {len(key_counts)} total keys "
            f"(min_count={min_count}, top_n={top_n}, min_frequency={min_frequency}%)"
        )

        if not filtered_keys:
            logger.warning("No keys passed the filter criteria")
            results = []
        else:
            # Perform alignment(s)
            results = self._align_filtered_keys(
                filtered_keys=filtered_keys,
                keys_to_barcodes=keys_to_barcodes,
                keys_to_sequences=keys_to_sequences,
                key_counts=key_counts,
                align_by_barcode=align_by_barcode,
                match=match,
                mismatch=mismatch,
                gap=gap,
                algorithm=algorithm,
            )

        # Write output
        self._write_aligned_output(
            results=results,
            output_file=output_file,
            key_col=key_col,
            barcode_col=barcode_col,
            align_by_barcode=align_by_barcode,
            sep=sep,
        )

        metrics = {
            "total_rows": total_rows,
            "total_keys": len(keys_to_barcodes),
            "filtered_keys": len(filtered_keys),
            "output_rows": len(results),
        }

        logger.info(
            f"Processed {len(keys_to_barcodes)} keys, "
            f"filtered to {len(filtered_keys)}, output {len(results)} rows"
        )
        return metrics

    def _filter_keys(
        self,
        key_counts: Dict[str, int],
        min_count: int,
        top_n: Optional[int],
        min_frequency: float,
    ) -> set:
        """Filter keys based on unique barcode count.

        Parameters
        ----------
        key_counts : Dict[str, int]
            Dictionary mapping keys to their unique barcode counts
        min_count : int
            Minimum unique barcode count threshold
        top_n : Optional[int]
            Keep only top N keys by unique barcode count
        min_frequency : float
            Minimum frequency percentage threshold

        Returns
        -------
        set
            Set of filtered keys
        """
        filtered_keys = set()
        
        # Calculate total unique barcode count across all keys
        total_unique_barcodes = sum(key_counts.values())
        
        # Apply filters
        for key, count in key_counts.items():
            if count < min_count:
                continue
            if min_frequency > 0 and total_unique_barcodes > 0:
                frequency = (count / total_unique_barcodes) * 100
                if frequency < min_frequency:
                    continue
            filtered_keys.add(key)
        
        # Apply top-n filter
        if top_n is not None and top_n > 0:
            sorted_keys = sorted(key_counts.items(), key=lambda x: x[1], reverse=True)
            top_keys = {key for key, _ in sorted_keys[:top_n]}
            filtered_keys &= top_keys
        
        return filtered_keys

    def _align_filtered_keys(
        self,
        filtered_keys: set,
        keys_to_barcodes: Dict[str, set],
        keys_to_sequences: Dict[str, List[str]],
        key_counts: Dict[str, int],
        align_by_barcode: bool,
        match: int,
        mismatch: int,
        gap: int,
        algorithm: int,
    ) -> List[Dict[str, Any]]:
        """Align filtered keys and return results.

        Parameters
        ----------
        filtered_keys : set
            Set of keys that passed filtering
        keys_to_barcodes : Dict[str, set]
            Dictionary mapping keys to sets of barcodes
        keys_to_sequences : Dict[str, List[str]]
            Dictionary mapping keys to lists of sequences
        key_counts : Dict[str, int]
            Dictionary mapping keys to their unique barcode counts
        align_by_barcode : bool
            If True, align keys separately grouped by barcode
        match : int
            Match score for alignment
        mismatch : int
            Mismatch penalty for alignment
        gap : int
            Gap penalty for alignment
        algorithm : int
            Alignment algorithm

        Returns
        -------
        List[Dict[str, Any]]
            List of result dictionaries with key, aligned_sequence, unique_barcode_count, and optionally barcode
        """
        results = []
        
        if align_by_barcode:
            # Group filtered keys by barcode
            # Keys with multiple barcodes will appear in multiple groups
            barcode_to_keys = defaultdict(list)
            for key in filtered_keys:
                for barcode in keys_to_barcodes[key]:
                    barcode_to_keys[barcode].append(key)
            
            # Align each barcode group separately
            for barcode, keys in tqdm(barcode_to_keys.items(), desc="Aligning by barcode"):
                # Collect sequences for this barcode group
                group_sequences = []
                key_to_sequence_indices = {}  # Track which sequences belong to which key
                
                for key in keys:
                    start_idx = len(group_sequences)
                    group_sequences.extend(keys_to_sequences[key])
                    end_idx = len(group_sequences)
                    key_to_sequence_indices[key] = (start_idx, end_idx)
                
                if not group_sequences:
                    continue
                
                # Align this group
                logger.debug(
                    f"Aligning {len(group_sequences)} sequences for barcode {barcode} "
                    f"({len(keys)} keys)"
                )
                aligned_sequences = align_umi_sequences(
                    sequences=group_sequences,
                    match=match,
                    mismatch=mismatch,
                    gap=gap,
                    algorithm=algorithm,
                )
                
                # Map aligned sequences back to keys
                # For simplicity, use the first aligned sequence for each key
                # (in practice, all sequences from the same key should align similarly)
                for key in keys:
                    start_idx, end_idx = key_to_sequence_indices[key]
                    # Use the first aligned sequence from this key's sequences
                    aligned_seq = aligned_sequences[start_idx] if aligned_sequences else ""
                    
                    results.append({
                        'key': key,
                        'barcode': barcode,
                        'aligned_sequence': aligned_seq,
                        'unique_barcode_count': key_counts[key]
                    })
        else:
            # Default mode: align all filtered keys together
            all_sequences = []
            key_to_sequence_indices = {}  # Track which sequences belong to which key
            
            for key in filtered_keys:
                start_idx = len(all_sequences)
                all_sequences.extend(keys_to_sequences[key])
                end_idx = len(all_sequences)
                key_to_sequence_indices[key] = (start_idx, end_idx)
            
            if not all_sequences:
                return results
            
            # Align all together
            logger.info(
                f"Aligning {len(all_sequences)} sequences from {len(filtered_keys)} keys "
                f"with match={match}, mismatch={mismatch}, gap={gap}, algorithm={algorithm}"
            )
            aligned_sequences = align_umi_sequences(
                sequences=all_sequences,
                match=match,
                mismatch=mismatch,
                gap=gap,
                algorithm=algorithm,
            )
            
            # Map aligned sequences back to keys
            for key in filtered_keys:
                start_idx, end_idx = key_to_sequence_indices[key]
                # Use the first aligned sequence from this key's sequences
                aligned_seq = aligned_sequences[start_idx] if aligned_sequences else ""
                
                results.append({
                    'key': key,
                    'aligned_sequence': aligned_seq,
                    'unique_barcode_count': key_counts[key]
                })
        
        return results

    def _write_aligned_output(
        self,
        results: List[Dict[str, Any]],
        output_file: Optional[str],
        key_col: str,
        barcode_col: str,
        align_by_barcode: bool,
        sep: str,
    ) -> None:
        """Write aligned sequences to CSV file or stdout.

        Parameters
        ----------
        results : List[Dict[str, Any]]
            List of result dictionaries with key, aligned_sequence, unique_barcode_count, and optionally barcode
        output_file : Optional[str]
            Path to output file (None for stdout)
        key_col : str
            Name of the key column
        barcode_col : str
            Name of the barcode column
        align_by_barcode : bool
            Whether barcode-grouped mode was used
        sep : str
            CSV separator character
        """
        if not results:
            logger.warning("No results to write")
            return
        
        if output_file:
            logger.info(f"Writing aligned sequences to {output_file}")
            # Create output directory if needed
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, "w", newline="") as f:
                writer = csv.writer(f, delimiter=sep)
                
                # Write header
                if align_by_barcode:
                    writer.writerow([key_col, barcode_col, "aligned_sequence", "unique_barcode_count"])
                else:
                    writer.writerow([key_col, "aligned_sequence", "unique_barcode_count"])
                
                # Write rows
                for result in results:
                    if align_by_barcode:
                        writer.writerow([
                            result['key'],
                            result['barcode'],
                            result['aligned_sequence'],
                            result['unique_barcode_count']
                        ])
                    else:
                        writer.writerow([
                            result['key'],
                            result['aligned_sequence'],
                            result['unique_barcode_count']
                        ])
        else:
            # Write to stdout
            logger.debug("Writing aligned sequences to stdout")
            
            # Write header
            if align_by_barcode:
                print(sep.join([key_col, barcode_col, "aligned_sequence", "unique_barcode_count"]))
            else:
                print(sep.join([key_col, "aligned_sequence", "unique_barcode_count"]))
            
            # Write rows
            for result in results:
                if align_by_barcode:
                    print(sep.join([
                        str(result['key']),
                        str(result['barcode']),
                        str(result['aligned_sequence']),
                        str(result['unique_barcode_count'])
                    ]))
                else:
                    print(sep.join([
                        str(result['key']),
                        str(result['aligned_sequence']),
                        str(result['unique_barcode_count'])
                    ]))

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
            "min_count": 0,
            "top_n": None,
            "min_frequency": 0.0,
            "align_by_barcode": False,
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.key_column:
            raise ValueError("Please provide --key-column")
        if not self.args.barcode_column:
            raise ValueError("Please provide --barcode-column")

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
                    key_col=self.args.key_column,
                    barcode_col=self.args.barcode_column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    match=self.args.match,
                    mismatch=self.args.mismatch,
                    gap=self.args.gap,
                    algorithm=self.args.algorithm,
                    min_count=self.args.min_count,
                    top_n=self.args.top_n,
                    min_frequency=self.args.min_frequency,
                    align_by_barcode=self.args.align_by_barcode,
                )

                logger.info(
                    f"Alignment complete. Processed {metrics.get('total_keys', 0)} keys, "
                    f"filtered to {metrics.get('filtered_keys', 0)}, "
                    f"output {metrics.get('output_rows', 0)} rows"
                )

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
                    key_col=self.args.key_column,
                    barcode_col=self.args.barcode_column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    match=self.args.match,
                    mismatch=self.args.mismatch,
                    gap=self.args.gap,
                    algorithm=self.args.algorithm,
                    min_count=self.args.min_count,
                    top_n=self.args.top_n,
                    min_frequency=self.args.min_frequency,
                    align_by_barcode=self.args.align_by_barcode,
                )

                logger.info(
                    f"Processed {os.path.basename(input_file)}: "
                    f"{metrics.get('total_keys', 0)} keys processed, "
                    f"{metrics.get('filtered_keys', 0)} filtered, "
                    f"{metrics.get('output_rows', 0)} rows output"
                )

            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise

        logger.info(
            f"Processing complete. Aligned sequences written to: {self.args.output_dir}"
        )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
