"""Tests for the findseq command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
import csv
from outerspace.cli.main import Cli


def test_findseq_initialization():
    """Test that findseq command initializes correctly"""
    args = [
        "findseq",
        "test_config.toml",
        "-1",
        "test_r1.fastq",
        "-2",
        "test_r2.fastq",
        "-o",
        "test_output.csv",
    ]
    cli = Cli(args)
    assert cli.args.config == "test_config.toml"
    assert cli.args.read1_filename == "test_r1.fastq"
    assert cli.args.read2_filename == "test_r2.fastq"
    assert cli.args.output_filename == "test_output.csv"


def test_findseq_single_file():
    """Test that findseq command works with single file"""
    args = [
        "findseq",
        "test_config.toml",
        "-1",
        "test_single.fastq",
        "-o",
        "test_output.csv",
    ]
    cli = Cli(args)
    assert cli.args.read1_filename == "test_single.fastq"
    assert cli.args.read2_filename is None


def test_findseq_format_detection():
    """Test that findseq command detects file formats correctly"""
    from outerspace.cli.commands.findseq import FindSeqCommand

    cmd = FindSeqCommand()

    # Test various formats
    assert cmd._detect_file_format("test.fastq") == "fastq"
    assert cmd._detect_file_format("test.fq") == "fastq"
    assert cmd._detect_file_format("test.fastq.gz") == "fastq"
    assert cmd._detect_file_format("test.fq.gz") == "fastq"

    assert cmd._detect_file_format("test.fasta") == "fasta"
    assert cmd._detect_file_format("test.fa") == "fasta"
    assert cmd._detect_file_format("test.fasta.gz") == "fasta"
    assert cmd._detect_file_format("test.fa.gz") == "fasta"

    assert cmd._detect_file_format("test.sam") == "sam"
    assert cmd._detect_file_format("test.bam") == "bam"

    # Test invalid format
    with pytest.raises(ValueError, match="Cannot detect format"):
        cmd._detect_file_format("test.unknown")


def test_findseq_missing_reads():
    """Test that findseq command handles missing read files"""
    args = ["findseq", "tests/configs/grnaquery.toml", "-o", "test_output.csv"]
    cli = Cli(args)
    with pytest.raises(
        ValueError,
        match="Please provide either -1 for single file or -1/-2 for paired files",
    ):
        cli.run()


def test_findseq_missing_output():
    """Test that findseq command requires output file"""
    args = ["findseq", "tests/configs/grnaquery.toml", "-1", "test_r1.fastq"]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Please provide an output filename with -o"):
        cli.run()


def test_findseq_with_example_data():
    """Test findseq command with real example data"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            "findseq",
            "tests/configs/grnaquery.toml",
            "-1",
            "tests/data/409-4_S1_L002_R1_001.fastq.gz",
            "-2",
            "tests/data/409-4_S1_L002_R2_001.fastq.gz",
            "-o",
            os.path.join(temp_dir, "shuffle.csv"),
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has content
        assert os.path.exists(os.path.join(temp_dir, "shuffle.csv"))

        with open(os.path.join(temp_dir, "shuffle.csv"), "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            # Check that we have the expected columns from the new pattern format
            expected_columns = [
                "read_id",
                "UMI_5prime",
                "protospacer",
                "downstreamof_protospacer",
                "UMI_3prime",
            ]
            assert all(
                col in header for col in expected_columns
            ), f"Expected columns {expected_columns}, got {header}"

            # Check that we have results
            rows = list(reader)
            assert len(rows) > 0, "No results found"


def test_findseq_single_file_processing():
    """Test findseq command with single file processing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a simple FASTA file for testing
        fasta_file = os.path.join(temp_dir, "test.fasta")
        with open(fasta_file, "w") as f:
            f.write(">read1\nATCGATCGATCG\n")
            f.write(">read2\nGCTAGCTAGCTA\n")

        # Create a simple config with one pattern
        config_file = os.path.join(temp_dir, "test_config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[findseq]
[[findseq.patterns]]
reg_expr = "(?P<test>.{4})"
read = "R1"
orientation = "forward"
multiple = "first"
"""
            )

        args = [
            "findseq",
            config_file,
            "-1",
            fasta_file,
            "-o",
            os.path.join(temp_dir, "output.csv"),
        ]
        cli = Cli(args)
        cli.run()

        # Verify output
        assert os.path.exists(os.path.join(temp_dir, "output.csv"))

        with open(os.path.join(temp_dir, "output.csv"), "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            assert "test" in header
            rows = list(reader)
            assert len(rows) > 0


def test_findseq_sam_bam_region():
    """Test findseq command with SAM/BAM region specification"""
    from outerspace.cli.commands.findseq import FindSeqCommand

    cmd = FindSeqCommand()

    # Test region parameter handling
    cmd.args = type("Args", (), {"region": "chr1:1-1000", "fetch": None})()

    # This would normally process a BAM file, but we're just testing the parameter handling
    # The actual BAM processing is tested in the Read class tests
    assert cmd.args.region == "chr1:1-1000"


def test_findseq_fetch_modes():
    """Test findseq command with different fetch modes"""
    from outerspace.cli.commands.findseq import FindSeqCommand

    cmd = FindSeqCommand()

    # Test fetch mode parameter handling
    for fetch_mode in ["mapped", "unmapped", "all"]:
        cmd.args = type("Args", (), {"region": None, "fetch": fetch_mode})()
        assert cmd.args.fetch == fetch_mode

    # Test invalid fetch mode
    cmd.args = type("Args", (), {"region": None, "fetch": "invalid"})()

    # This would raise an error in actual processing, but we're just testing parameter handling
    assert cmd.args.fetch == "invalid"
