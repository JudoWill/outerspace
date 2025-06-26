"""Sequence extraction command for pattern matching.

This module provides the FindSeqCommand class for extracting sequences from various
file formats based on user-defined patterns. It supports FASTQ, FASTA, SAM, and BAM
formats with both single-end and paired-end processing.
"""

import csv
import logging
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Union

from outerspace.cli.commands.base import BaseCommand
from outerspace.config import Cfg
from outerspace.pattern import Pattern
from outerspace.read import Read

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def run(
    file_pathname_cfg: str,
    file_pathname_fastq1: str,
    file_pathname_fastq2: Optional[str] = None,
    filename_csv: Optional[str] = None,
) -> None:
    """Parse reads, find user-given patterns, save sequence matches to CSV.

    This is a convenience wrapper around the findseq command for programmatic use.

    Parameters
    ----------
    file_pathname_cfg : str
        Path to configuration file with search patterns
    file_pathname_fastq1 : str
        Path to read 1 file (FASTQ, FASTA, SAM, BAM)
    file_pathname_fastq2 : Optional[str], default=None
        Path to read 2 file for paired-end data
    filename_csv : Optional[str], default=None
        Path to output CSV file for results

    Raises
    ------
    ValueError
        If any input files don't exist
    """
    if not Path(file_pathname_cfg).exists():
        raise ValueError(f"Configuration file not found: {file_pathname_cfg}")
    if not Path(file_pathname_fastq1).exists():
        raise ValueError(f"Read 1 file not found: {file_pathname_fastq1}")
    if file_pathname_fastq2 and not Path(file_pathname_fastq2).exists():
        raise ValueError(f"Read 2 file not found: {file_pathname_fastq2}")

    # Create command object
    cmd = FindSeqCommand()

    # Set up arguments
    cmd.args = type(
        "Args",
        (),
        {
            "config": file_pathname_cfg,
            "read1_filename": file_pathname_fastq1,
            "read2_filename": file_pathname_fastq2,
            "output_filename": filename_csv,
        },
    )()

    # Run the command
    cmd.run()


class FindSeqCommand(BaseCommand):
    """Command for extracting sequences from various file formats.

    This command processes sequence files (FASTQ, FASTA, SAM, BAM) and extracts
    sequences that match user-defined patterns. It supports both single-end and
    paired-end data with comprehensive pattern matching capabilities.
    """

    def _init_parser(self, subparsers) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparsers
            Subparser group to add command arguments to
        """
        parser = subparsers.add_parser(
            "findseq",
            help="Extract sequences from files based on configuration patterns",
        )
        parser.add_argument("config", help="Configuration file with search patterns")
        parser.add_argument(
            "-1",
            "--read1_filename",
            help="Input file for read 1 (FASTQ, FASTA, SAM, BAM) or single read file",
        )
        parser.add_argument(
            "-2",
            "--read2_filename",
            help="Input file for read 2 (FASTQ, FASTA, SAM, BAM) for paired reads",
        )
        parser.add_argument("-o", "--output_filename", help="Output CSV file name")
        parser.add_argument(
            "--region", help='SAM/BAM region specification (e.g., "chr1:1-1000")'
        )
        parser.add_argument(
            "--fetch", help="SAM/BAM fetch mode (mapped, unmapped, all)"
        )

    def _detect_file_format(self, filename: str) -> str:
        """Auto-detect file format based on extension.

        Parameters
        ----------
        filename : str
            Path to the file to analyze

        Returns
        -------
        str
            Detected format ('fastq', 'fasta', 'sam', 'bam')

        Raises
        ------
        ValueError
            If file format cannot be detected
        """
        key = filename.lower().replace(".gz", "")
        if key.endswith((".fastq", ".fq")):
            return "fastq"
        elif key.endswith((".fasta", ".fa")):
            return "fasta"
        elif key.endswith(".sam"):
            return "sam"
        elif key.endswith(".bam"):
            return "bam"
        else:
            raise ValueError(
                f"Cannot detect format for file: {filename}. "
                f"Please use a file with .fastq, .fasta, .sam, or .bam extension."
            )

    def _process_single_file(self, filename: str) -> Generator[Read, None, None]:
        """Process a single file of any supported format.

        Parameters
        ----------
        filename : str
            Path to the file to process

        Yields
        ------
        Read
            Read objects from the file

        Raises
        ------
        ValueError
            If file format is unsupported or fetch mode is invalid
        """
        format_type = self._detect_file_format(filename)
        logger.info(f"Processing single file: {filename} (format: {format_type})")

        if format_type in ["sam", "bam"]:
            # Handle SAM/BAM with optional region and fetch parameters
            fetch_param = None
            if self.args.region:
                fetch_param = self.args.region
                logger.debug(f"Using region: {self.args.region}")
            elif self.args.fetch:
                if self.args.fetch == "mapped":
                    fetch_param = True
                elif self.args.fetch == "unmapped":
                    fetch_param = False
                elif self.args.fetch == "all":
                    fetch_param = None
                else:
                    raise ValueError(f"Invalid fetch mode: {self.args.fetch}")
                logger.debug(f"Using fetch mode: {self.args.fetch}")

            yield from Read.from_bam(filename, fetch=fetch_param)
        elif format_type == "fasta":
            yield from Read.from_fasta(filename)
        elif format_type == "fastq":
            yield from Read.from_fastq(filename)
        else:
            raise ValueError(f"Unsupported format: {format_type}")

    def _process_paired_files(
        self, file1: str, file2: str
    ) -> Generator[Read, None, None]:
        """Process paired files of any supported format.

        Parameters
        ----------
        file1 : str
            Path to first file in pair
        file2 : str
            Path to second file in pair

        Yields
        ------
        Read
            Read objects from the paired files

        Raises
        ------
        ValueError
            If paired files have different formats
        """
        format1 = self._detect_file_format(file1)
        format2 = self._detect_file_format(file2)

        if format1 != format2:
            raise ValueError(
                f"Paired files must have the same format. Got {format1} and {format2}"
            )

        logger.info(f"Processing paired files: {file1}, {file2} (format: {format1})")
        yield from Read.from_paired_fastx(file1, file2, format=format1)

    def _process_reads_with_patterns(
        self, reads: Generator[Read, None, None]
    ) -> List[Dict[str, Any]]:
        """Process reads using Pattern objects and return results.

        Parameters
        ----------
        reads : Generator[Read, None, None]
            Generator of Read objects to process

        Returns
        -------
        List[Dict[str, Any]]
            List of dictionaries containing pattern match results

        Notes
        -----
        This method applies all configured patterns to each read and collects
        captured groups from successful matches. It handles the 'multiple'
        parameter for each pattern (first, last, all).
        """
        patterns = self._parse_patterns_from_config()
        logger.info(f"Processing reads with {len(patterns)} patterns")

        results = []
        read_count = 0
        match_count = 0

        for read in reads:
            read_count += 1
            read_results = {}

            for pattern in patterns:
                hits = pattern.search(read)
                if hits:
                    # Handle multiple, first, last logic
                    if pattern.multiple == "first":
                        hits = [hits] if not isinstance(hits, list) else [hits[0]]
                    elif pattern.multiple == "last":
                        hits = [hits] if not isinstance(hits, list) else [hits[-1]]
                    # hits is already a list for 'all'

                    for hit in hits:
                        read_results.update(hit.captured)

            if read_results:
                read_results["read_id"] = read.name
                results.append(read_results)
                match_count += 1

        logger.info(f"Processed {read_count} reads, found {match_count} matches")
        return results

    def _write_results_to_csv(
        self, results: List[Dict[str, Any]], output_filename: str
    ) -> None:
        """Write results to CSV file.

        Parameters
        ----------
        results : List[Dict[str, Any]]
            List of result dictionaries to write
        output_filename : str
            Path to output CSV file

        Notes
        -----
        This method automatically determines all unique keys from the results
        and creates a CSV with consistent column ordering.
        """
        if not results:
            logger.warning("No matches found")
            return

        # Get all unique keys from all results
        all_keys = set()
        for result in results:
            all_keys.update(result.keys())

        # Sort keys for consistent output
        fieldnames = sorted(all_keys)

        logger.info(f"Writing {len(results)} results to {output_filename}")

        with open(output_filename, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        logger.info(f"Results written to {output_filename}")

    def run(self) -> None:
        """Run the findseq command.

        This method orchestrates the sequence extraction process, handling
        configuration loading, file processing, pattern matching, and result
        output with comprehensive error handling.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        logger.info("Starting sequence extraction process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "read1_filename": None,
            "read2_filename": None,
            "output_filename": None,
            "region": None,
            "fetch": None,
        }
        self._merge_config_and_args(defaults)

        if self.args.config is None:
            raise ValueError("Please provide a config filename")

        # Check for required output file
        if not self.args.output_filename:
            raise ValueError("Please provide an output filename with -o")

        # Process based on input files
        if self.args.read1_filename and self.args.read2_filename:
            # Paired files
            self._chk_exists([self.args.read1_filename, self.args.read2_filename])
            reads = self._process_paired_files(
                self.args.read1_filename, self.args.read2_filename
            )
        elif self.args.read1_filename:
            # Single file
            self._chk_exists([self.args.read1_filename])
            reads = self._process_single_file(self.args.read1_filename)
        else:
            raise ValueError(
                "Please provide either -1 for single file or -1/-2 for paired files"
            )

        # Process reads with patterns
        results = self._process_reads_with_patterns(reads)

        # Write results
        self._write_results_to_csv(results, self.args.output_filename)

        logger.info("Sequence extraction process completed")


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
