"""Sequence extraction command for pattern matching.

This module provides the FindSeqCommand class for extracting sequences from various
file formats based on user-defined patterns. It supports FASTQ, FASTA, SAM, and BAM
formats with both single-end and paired-end processing.
"""

import csv
import logging
import os
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from itertools import islice
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Union

from outerspace.cli.commands.base import BaseCommand
from outerspace.config import Cfg
from outerspace.pattern import Hit, Pattern
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
    long_format: bool = False,
    matches_only: bool = False,
    skip_unmapped: bool = False,
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
    long_format : bool, default=False
        Output in long format (one row per pattern match)
    matches_only : bool, default=False
        Only output reads that have at least one pattern match
    skip_unmapped : bool, default=False
        Skip unmapped reads in SAM/BAM files (based on SAM flag 0x4)

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
            "long_format": long_format,
            "matches_only": matches_only,
            "skip_unmapped": skip_unmapped,
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

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "findseq",
            help="Extract sequences from files based on configuration patterns",
        )
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
        parser.add_argument(
            "--long-format",
            action="store_true",
            help="Output in long format (one row per pattern match instead of one row per read)",
        )
        parser.add_argument(
            "--matches-only",
            action="store_true",
            help="Only output reads that have at least one pattern match",
        )
        parser.add_argument(
            "--threads",
            type=int,
            default=None,
            help="Number of threads for parallel processing (default: auto-detect)",
        )
        parser.add_argument(
            "--skip-unmapped",
            action="store_true",
            help="Skip unmapped reads in SAM/BAM files (based on SAM flag 0x4)",
        )
        parser.add_argument(
            "--max-reads",
            type=int,
            default=None,
            help="Maximum number of reads to process (default: all)",
        )
        self._add_common_args(parser)

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

            yield from Read.from_bam(filename, fetch=fetch_param, skip_unmapped=self.args.skip_unmapped)
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

    def _chunk_reads(
        self, reads: Generator[Read, None, None], chunk_size: int = 100
    ) -> Generator[List[Read], None, None]:
        """Convert read generator into chunks for parallel processing.

        Parameters
        ----------
        reads : Generator[Read, None, None]
            Generator of Read objects
        chunk_size : int, default=1000
            Number of reads per chunk

        Yields
        ------
        List[Read]
            Lists of Read objects (chunks)
        """
        iterator = iter(reads)
        while True:
            chunk = list(islice(iterator, chunk_size))
            if not chunk:
                break
            yield chunk

    def _process_read_chunk(
        self, read_chunk: List[Read], patterns: List[Pattern]
    ) -> List[tuple[str, Any]]:
        """Process a chunk of reads with all patterns.

        Parameters
        ----------
        read_chunk : List[Read]
            List of Read objects to process
        patterns : List[Pattern]
            List of Pattern objects to apply

        Returns
        -------
        List[tuple[str, Any]]
            List of (readname, hit) tuples for all matches in the chunk
        """
        results = []
        
        for read in read_chunk:
            read_has_matches = False

            for pattern in patterns:
                hits = pattern.search(read)
                if hits:
                    read_has_matches = True
                    # Handle multiple, first, last logic
                    if pattern.multiple == "first":
                        hits = [hits] if not isinstance(hits, list) else [hits[0]]
                    elif pattern.multiple == "last":
                        hits = [hits] if not isinstance(hits, list) else [hits[-1]]
                    # hits is already a list for 'all'

                    for hit in hits:
                        results.append((read.name, hit))

            # If not matches-only and no matches found, yield empty hit for wide format
            if not read_has_matches and not self.args.matches_only:
                # Create an empty hit object for reads with no matches
                empty_hit = Hit(
                    start=-1,
                    end=-1,
                    match="",
                    orientation="",
                    captured={},
                    pattern=None,
                )
                results.append((read.name, empty_hit))
        
        return results

    def _process_reads_with_patterns(
        self, reads: Generator[Read, None, None]
    ) -> Generator[tuple[str, Any], None, None]:
        """Process reads using Pattern objects and yield (readname, hit) tuples.

        This method uses multi-threading to process reads in parallel for improved performance.
        Progress tracking is handled by wrapping the input reads generator.

        Parameters
        ----------
        reads : Generator[Read, None, None]
            Generator of Read objects to process

        Yields
        ------
        tuple[str, Any]
            Tuples of (readname, hit) where hit contains pattern match information

        Notes
        -----
        This method applies all configured patterns to each read and yields
        (readname, hit) tuples for each pattern match. It handles the 'multiple'
        parameter for each pattern (first, last, all).
        """
        patterns = self._parse_patterns_from_config()
        logger.info(f"Processing reads with {len(patterns)} patterns")

        # Determine number of threads
        max_workers = self.args.threads
        if max_workers is None:
            max_workers = min(os.cpu_count() or 1, 8)  # Default: auto-detect, max 8
        
        logger.info(f"Using {max_workers} threads for processing")

        # If only 1 thread, use single-threaded approach for simplicity
        if max_workers == 1:
            yield from self._process_reads_single_threaded(reads, patterns)
            return

        # Multi-threaded processing
        read_count = 0
        chunk_size = 1000  # Configurable if needed later

        # Convert reads to chunks
        read_chunks = self._chunk_reads(reads, chunk_size)

        # Process chunks in parallel using executor.map to maintain order
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Create a partial function with patterns bound
            process_chunk_with_patterns = partial(self._process_read_chunk, patterns=patterns)
            
            # Process chunks and yield results in order
            for chunk_results in executor.map(process_chunk_with_patterns, read_chunks):
                read_count += len([r for r in chunk_results if r[1].pattern is not None or not self.args.matches_only])
                for result in chunk_results:
                    yield result

        logger.info(f"Processed approximately {read_count} reads")

    def _process_reads_single_threaded(
        self, reads: Generator[Read, None, None], patterns: List[Pattern]
    ) -> Generator[tuple[str, Any], None, None]:
        """Process reads in single-threaded mode (original implementation).

        Parameters
        ----------
        reads : Generator[Read, None, None]
            Generator of Read objects to process
        patterns : List[Pattern]
            List of Pattern objects to apply

        Yields
        ------
        tuple[str, Any]
            Tuples of (readname, hit) where hit contains pattern match information
        """
        read_count = 0
        match_count = 0

        for read in reads:
            read_count += 1
            read_has_matches = False

            for pattern in patterns:
                hits = pattern.search(read)
                if hits:
                    read_has_matches = True
                    # Handle multiple, first, last logic
                    if pattern.multiple == "first":
                        hits = [hits] if not isinstance(hits, list) else [hits[0]]
                    elif pattern.multiple == "last":
                        hits = [hits] if not isinstance(hits, list) else [hits[-1]]
                    # hits is already a list for 'all'

                    for hit in hits:
                        yield (read.name, hit)
                        match_count += 1

            # If not matches-only and no matches found, yield empty hit for wide format
            if not read_has_matches and not self.args.matches_only:
                # Create an empty hit object for reads with no matches
                empty_hit = Hit(
                    start=-1,
                    end=-1,
                    match="",
                    orientation="",
                    captured={},
                    pattern=None,
                )
                yield (read.name, empty_hit)

        logger.info(f"Processed {read_count} reads, found {match_count} matches")

    def _write_results_to_csv_wide(
        self, results: Generator[tuple[str, Any], None, None], output_filename: str
    ) -> None:
        """Write results to CSV file in wide format (one row per read).

        Parameters
        ----------
        results : Generator[tuple[str, Any], None, None]
            Generator of (readname, hit) tuples
        output_filename : str
            Path to output CSV file

        Notes
        -----
        This method combines all pattern matches for each read into a single row.
        """
        # Collect all results to group by read_id
        read_data = {}
        all_keys = set()

        for readname, hit in results:
            if readname not in read_data:
                read_data[readname] = {"read_id": readname}

            # Add captured groups to the read's data
            for key, value in hit.captured.items():
                read_data[readname][key] = value
                all_keys.add(key)

        if not read_data:
            logger.warning("No matches found")
            return

        # Sort keys for consistent output, but ensure read_id comes first
        fieldnames = ["read_id"] + sorted(all_keys)

        logger.info(f"Writing {len(read_data)} results to {output_filename}")

        with open(output_filename, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(read_data.values())

        logger.info(f"Results written to {output_filename}")

    def _write_results_to_csv_long(
        self, results: Generator[tuple[str, Any], None, None], output_filename: str
    ) -> None:
        """Write results to CSV file in long format (one row per pattern match).

        Parameters
        ----------
        results : Generator[tuple[str, Any], None, None]
            Generator of (readname, hit) tuples
        output_filename : str
            Path to output CSV file

        Notes
        -----
        This method creates one row per pattern match with columns: read_id, pattern_name, match, start.
        """
        fieldnames = ["read_id", "pattern_name", "match", "start"]

        logger.info(f"Writing long format results to {output_filename}")

        with open(output_filename, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            row_count = 0
            for readname, hit in results:
                # Skip empty hits (reads with no matches)
                if hit.pattern is None:
                    continue
                if hit.captured:
                    # If there are captured groups, write one row per group
                    for key, value in hit.captured.items():
                        row = {
                            "read_id": readname,
                            "pattern_name": key,
                            "match": value,
                            "start": getattr(hit, "start", -1),
                        }
                        writer.writerow(row)
                        row_count += 1
                else:
                    # If no captured groups, write one row per pattern match
                    row = {
                        "read_id": readname,
                        "pattern_name": hit.pattern.name,
                        "match": "",
                        "start": -1,
                    }
                    writer.writerow(row)
                    row_count += 1

        logger.info(f"Results written to {output_filename} ({row_count} rows)")

    def _write_results_to_csv(
        self, results: Generator[tuple[str, Any], None, None], output_filename: str
    ) -> None:
        """Write results to CSV file in the appropriate format.

        Parameters
        ----------
        results : Generator[tuple[str, Any], None, None]
            Generator of (readname, hit) tuples
        output_filename : str
            Path to output CSV file
        """
        if self.args.long_format:
            self._write_results_to_csv_long(results, output_filename)
        else:
            self._write_results_to_csv_wide(results, output_filename)

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
            "long_format": False,
            "matches_only": False,
            "threads": None,
            "skip_unmapped": False,
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

        # Wrap reads with progress bar for cleaner progress tracking
        reads_with_progress = self._progbar_iterable(reads, desc="Processing reads")

        # Process reads with patterns (multi-threaded)
        results = self._process_reads_with_patterns(reads_with_progress)

        # Write results
        self._write_results_to_csv(results, self.args.output_filename)

        logger.info("Sequence extraction process completed")


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
