"""Read handling functionality for sequence data files.

This module provides the Read class for representing sequence reads and methods
for reading from various file formats including FASTA, FASTQ, and BAM/SAM files.
It supports both single-end and paired-end sequencing data.
"""

import logging
from typing import Generator, Optional
from itertools import zip_longest

import pyfastx
import pysam

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class Read:
    """Container for sequence read information.

    This class represents a single sequence read with its sequence data,
    read pair designation, and optional name. It provides methods for
    reverse complement calculation and various file format readers.
    """

    # Valid read pair designations
    VALID_PAIRS = {"R1", "R2"}

    def __init__(self, seq: str, pair: str, name: Optional[str] = None) -> None:
        """Initialize a Read object.

        Parameters
        ----------
        seq : str
            The DNA/RNA sequence string
        pair : str
            Read pair designation ('R1' or 'R2')
        name : Optional[str], default=None
            Optional name/identifier for the read

        Raises
        ------
        ValueError
            If pair is not a valid read pair designation
        """
        if pair not in self.VALID_PAIRS:
            raise ValueError(f"Invalid pair: {pair}. Must be one of {self.VALID_PAIRS}")

        self.name = name
        self.seq = seq
        self.pair = pair

        logger.debug(f"Initialized read: {name} ({pair}) with {len(seq)} bp")

    @property
    def seq_rc(self) -> str:
        """Get the reverse complement of the sequence.

        Returns
        -------
        str
            Reverse complement sequence

        Notes
        -----
        This property computes the reverse complement by reversing the sequence
        and translating A↔T and C↔G (case-insensitive).
        """
        return self.seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))

    def __str__(self) -> str:
        """Return string representation of the read.

        Returns
        -------
        str
            Formatted string showing sequence and pair designation
        """
        return f"{self.seq}\n{self.pair}"

    def __repr__(self) -> str:
        """Return detailed string representation for debugging.

        Returns
        -------
        str
            Same as __str__ for this class
        """
        return f"{self.seq}\n{self.pair}"

    @classmethod
    def from_fasta(
        cls, fasta_file: str, read_pair: str = "R1"
    ) -> Generator["Read", None, None]:
        """Iterate Read objects from a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to the FASTA file
        read_pair : str, default='R1'
            Read pair designation for all reads in the file

        Yields
        ------
        Read
            Read objects for each sequence in the FASTA file

        Notes
        -----
        Uses pyfastx for efficient FASTA parsing without building an index.
        """
        try:
            for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
                yield cls(seq=str(seq), name=name, pair=read_pair)
        except Exception as e:
            logger.error(f"Error reading FASTA file {fasta_file}: {e}")
            raise

    @classmethod
    def from_fastq(
        cls, fastq_file: str, read_pair: str = "R1"
    ) -> Generator["Read", None, None]:
        """Iterate Read objects from a FASTQ file.

        Parameters
        ----------
        fastq_file : str
            Path to the FASTQ file
        read_pair : str, default='R1'
            Read pair designation for all reads in the file

        Yields
        ------
        Read
            Read objects for each sequence in the FASTQ file

        Notes
        -----
        Uses pyfastx for efficient FASTQ parsing without building an index.
        Quality scores are ignored in the current implementation.
        """
        try:
            for name, seq, _ in pyfastx.Fastq(fastq_file, build_index=False):
                yield cls(seq=str(seq), name=name, pair=read_pair)
        except Exception as e:
            logger.error(f"Error reading FASTQ file {fastq_file}: {e}")
            raise

    @classmethod
    def from_paired_fastx(
        cls, fastx_file1: str, fastx_file2: str, format: str = "fasta"
    ) -> Generator["Read", None, None]:
        """Iterate Read objects from paired FASTX files.

        This method reads from two files and yields reads in alternating order
        (R1, R2, R1, R2, ...) to maintain pairing information.

        Parameters
        ----------
        fastx_file1 : str
            Path to the first FASTX file (R1 reads)
        fastx_file2 : str
            Path to the second FASTX file (R2 reads)
        format : str, default='fasta'
            File format ('fasta' or 'fastq')

        Yields
        ------
        Read
            Read objects alternating between R1 and R2 reads

        Raises
        ------
        ValueError
            If format is not supported

        Notes
        -----
        Uses zip_longest to handle cases where files have different numbers
        of reads. Missing reads are skipped.
        """
        format_lower = format.lower()

        if format_lower == "fasta":
            fx1 = cls.from_fasta(fastx_file1, read_pair="R1")
            fx2 = cls.from_fasta(fastx_file2, read_pair="R2")
        elif format_lower == "fastq":
            fx1 = cls.from_fastq(fastx_file1, read_pair="R1")
            fx2 = cls.from_fastq(fastx_file2, read_pair="R2")
        else:
            raise ValueError(
                f"Unsupported format: {format}. Must be 'fasta' or 'fastq'"
            )

        logger.info(f"Reading paired {format} files: {fastx_file1}, {fastx_file2}")

        # Iterate through both files in alternating fashion
        for read1, read2 in zip_longest(fx1, fx2):
            if read1 is not None:
                yield read1
            if read2 is not None:
                yield read2

    @classmethod
    def from_bam(
        cls, bam_file: str, fetch: Optional[str] = None
    ) -> Generator["Read", None, None]:
        """Iterate Read objects from a SAM/BAM file.

        Parameters
        ----------
        bam_file : str
            Path to the SAM/BAM file
        fetch : Optional[str], default=None
            Optional region specification for targeted reading.
            Format: "chr1:1-1000" or "chr1"

        Yields
        ------
        Read
            Read objects for each alignment in the file

        Notes
        -----
        Automatically detects file format based on extension (.sam vs .bam).
        Read pair designation is inferred from SAM flags:
        - Bit 6 (64): First in pair (R1)
        - Bit 7 (128): Second in pair (R2)
        - Neither: Single-end read (assigned as R1)
        """
        try:
            # Determine file format and open accordingly
            if bam_file.endswith(".sam"):
                mode = "r"  # SAM file
                logger.debug(f"Opening SAM file: {bam_file}")
            else:
                mode = "rb"  # BAM file
                logger.debug(f"Opening BAM file: {bam_file}")

            with pysam.AlignmentFile(bam_file, mode) as bam:
                if fetch is not None:
                    reads = cls._fetch_region(bam, fetch)
                else:
                    reads = bam

                for read in reads:
                    # Infer read pair from SAM flag
                    read_pair = cls._infer_read_pair(read)

                    yield cls(
                        seq=read.query_sequence, name=read.query_name, pair=read_pair
                    )

        except Exception as e:
            logger.error(f"Error reading BAM/SAM file {bam_file}: {e}")
            raise

    @staticmethod
    def _fetch_region(bam: pysam.AlignmentFile, fetch: str) -> pysam.IteratorRow:
        """Fetch reads from a specific genomic region.

        Parameters
        ----------
        bam : pysam.AlignmentFile
            Open BAM/SAM file object
        fetch : str
            Region specification in format "chr1:1-1000" or "chr1"

        Returns
        -------
        pysam.IteratorRowRegion
            Iterator over reads in the specified region

        Notes
        -----
        Coordinates are converted from 1-based (user) to 0-based (pysam) as needed.
        """
        if ":" in fetch:
            # Region with coordinates: "chr1:1-10"
            contig, coords = fetch.split(":", 1)
            if "-" in coords:
                start_str, stop_str = coords.split("-", 1)
                start = int(start_str) - 1  # Convert to 0-based
                stop = int(stop_str)  # Keep as 1-based for pysam
            else:
                start = int(coords) - 1  # Single position
                stop = start + 1
        else:
            # Just contig name: "chr1"
            contig = fetch
            start = None
            stop = None

        logger.debug(f"Fetching region: {contig}:{start}-{stop}")
        return bam.fetch(contig, start, stop)

    @staticmethod
    def _infer_read_pair(read: pysam.AlignedSegment) -> str:
        """Infer read pair designation from SAM flags.

        Parameters
        ----------
        read : pysam.AlignedSegment
            Aligned read object

        Returns
        -------
        str
            Read pair designation ('R1' or 'R2')

        Notes
        -----
        Uses SAM flags to determine pairing:
        - Bit 6 (64): First in pair
        - Bit 7 (128): Second in pair
        - Neither: Single-end read (assigned as R1)
        """
        if read.flag & 64:  # First in pair
            return "R1"
        elif read.flag & 128:  # Second in pair
            return "R2"
        else:  # Single end or unpaired
            return "R1"


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
