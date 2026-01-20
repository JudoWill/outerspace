"""Tests for the align command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli


def test_align_initialization():
    """Test that align command initializes correctly"""
    args = [
        "align",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--column",
        "umi",
    ]
    cli = Cli(args)
    assert cli.args.input_dir == "test_input"
    assert cli.args.output_dir == "test_output"
    assert cli.args.column == "umi"
    assert cli.args.sep == ","
    assert cli.args.match == 5
    assert cli.args.mismatch == -4
    assert cli.args.gap == -8
    assert cli.args.algorithm == 1


def test_align_single_file_initialization():
    """Test that align command initializes correctly in single file mode"""
    args = [
        "align",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.txt",
        "--column",
        "umi",
    ]
    cli = Cli(args)
    assert cli.args.input_file == "test_input.csv"
    assert cli.args.output_file == "test_output.txt"
    assert cli.args.column == "umi"
    assert cli.args.sep == ","
    assert cli.args.match == 5
    assert cli.args.mismatch == -4
    assert cli.args.gap == -8
    assert cli.args.algorithm == 1


def test_align_missing_input():
    """Test that align command handles missing input"""
    args = [
        "align",
        "--output-dir",
        "test_output",
        "--column",
        "umi",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_align_missing_column():
    """Test that align command handles missing column"""
    args = [
        "align",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.txt",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_align_missing_output_dir_for_directory_mode():
    """Test that align command requires output-dir for directory mode"""
    args = [
        "align",
        "--input-dir",
        "test_input",
        "--column",
        "umi",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Please provide --output-dir"):
        cli.run()


def test_align_single_file_nonexistent_input():
    """Test that align command handles nonexistent input file"""
    args = [
        "align",
        "--input-file",
        "nonexistent.csv",
        "--output-file",
        "test_output.txt",
        "--column",
        "umi",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Input file not found"):
        cli.run()


def test_align_directory_nonexistent_input():
    """Test that align command handles nonexistent input directory"""
    args = [
        "align",
        "--input-dir",
        "nonexistent_dir",
        "--output-dir",
        "test_output",
        "--column",
        "umi",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Input directory not found"):
        cli.run()


def test_align_empty_input_directory():
    """Test that align command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            "align",
            "--input-dir",
            temp_dir,
            "--output-dir",
            os.path.join(temp_dir, "output"),
            "--column",
            "umi",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="No CSV files found"):
            cli.run()


def test_align_invalid_column():
    """Test that align command handles invalid column name"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("header1,header2\nvalue1,value2\n")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            os.path.join(temp_dir, "output.txt"),
            "--column",
            "nonexistent_column",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Column 'nonexistent_column' not found"):
            cli.run()


def test_align_single_file_basic():
    """Test basic single file alignment functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with UMI sequences
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")
            f.write("AACTATA,data3\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)

        # Read and verify aligned sequences
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]

        # Should have 3 aligned sequences
        assert len(aligned_lines) == 3

        # All sequences should have the same length (aligned)
        aligned_lengths = set(len(seq) for seq in aligned_lines)
        assert len(aligned_lengths) == 1, "All aligned sequences should have the same length"


def test_align_single_file_stdout():
    """Test single file alignment output to stdout"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")

        args = [
            "align",
            "--input-file",
            input_file,
            "--column",
            "umi",
        ]
        cli = Cli(args)
        # Should not raise an error (outputs to stdout)
        cli.run()


def test_align_single_file_with_parameters():
    """Test single file alignment with custom parameters"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")
            f.write("AACTATA,data3\n")

        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
            "--match",
            "10",
            "--mismatch",
            "-2",
            "--gap",
            "-5",
            "--algorithm",
            "0",  # local alignment
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]
        assert len(aligned_lines) == 3


def test_align_directory_mode():
    """Test directory batch processing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create input directory with CSV files
        input_dir = os.path.join(temp_dir, "input")
        os.makedirs(input_dir)

        # Create first CSV file
        file1 = os.path.join(input_dir, "file1.csv")
        with open(file1, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")

        # Create second CSV file
        file2 = os.path.join(input_dir, "file2.csv")
        with open(file2, "w") as f:
            f.write("umi,other\n")
            f.write("GGGTTTAA,data3\n")
            f.write("GGGTTTAG,data4\n")

        # Create output directory
        output_dir = os.path.join(temp_dir, "output")

        args = [
            "align",
            "--input-dir",
            input_dir,
            "--output-dir",
            output_dir,
            "--column",
            "umi",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output files exist
        output_file1 = os.path.join(output_dir, "file1.csv")
        output_file2 = os.path.join(output_dir, "file2.csv")

        assert os.path.exists(output_file1)
        assert os.path.exists(output_file2)

        # Verify content
        with open(output_file1, "r") as f:
            lines1 = [line.strip() for line in f.readlines() if line.strip()]
            assert len(lines1) == 2

        with open(output_file2, "r") as f:
            lines2 = [line.strip() for line in f.readlines() if line.strip()]
            assert len(lines2) == 2


def test_align_empty_sequences():
    """Test alignment with empty sequences in CSV"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with some empty sequences
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write(",data2\n")  # Empty UMI
            f.write("AACTATA,data3\n")
            f.write(",data4\n")  # Empty UMI

        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
        ]
        cli = Cli(args)
        cli.run()

        # Should only align non-empty sequences
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]
        # Should have 2 aligned sequences (empty ones filtered out)
        assert len(aligned_lines) == 2


def test_align_single_sequence():
    """Test alignment with only one sequence"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with single sequence
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")

        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
        ]
        cli = Cli(args)
        cli.run()

        # Single sequence should be returned as-is
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]
        assert len(aligned_lines) == 1
        assert aligned_lines[0] == "AACTTATA"


def test_align_row_limit():
    """Test alignment with row limit"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")
            f.write("AACTATA,data3\n")
            f.write("GGGTTTAA,data4\n")

        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
            "--row-limit",
            "3",  # Only process first 3 data rows (3 sequences)
        ]
        cli = Cli(args)
        cli.run()

        # Should only align sequences from first 3 data rows
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]
        assert len(aligned_lines) == 3


def test_align_different_algorithms():
    """Test alignment with different algorithm types"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")
            f.write("AACTTATG,data2\n")
            f.write("AACTATA,data3\n")

        # Test local algorithm (0)
        output_file1 = os.path.join(temp_dir, "output_local.txt")
        args1 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file1,
            "--column",
            "umi",
            "--algorithm",
            "0",
        ]
        cli1 = Cli(args1)
        cli1.run()

        # Test global algorithm (1)
        output_file2 = os.path.join(temp_dir, "output_global.txt")
        args2 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file2,
            "--column",
            "umi",
            "--algorithm",
            "1",
        ]
        cli2 = Cli(args2)
        cli2.run()

        # Test semi-global algorithm (2)
        output_file3 = os.path.join(temp_dir, "output_semiglobal.txt")
        args3 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file3,
            "--column",
            "umi",
            "--algorithm",
            "2",
        ]
        cli3 = Cli(args3)
        cli3.run()

        # All should produce valid output
        assert os.path.exists(output_file1)
        assert os.path.exists(output_file2)
        assert os.path.exists(output_file3)


def test_align_invalid_output_dir_with_single_file():
    """Test that using --output-dir with --input-file raises error"""
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi,other\n")
            f.write("AACTTATA,data1\n")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-dir",
            os.path.join(temp_dir, "output"),
            "--column",
            "umi",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Cannot use --output-dir with --input-file"):
            cli.run()


def test_align_custom_separator():
    """Test alignment with custom CSV separator"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with tab separator
        input_file = os.path.join(temp_dir, "test.tsv")
        with open(input_file, "w") as f:
            f.write("umi\tother\n")
            f.write("AACTTATA\tdata1\n")
            f.write("AACTTATG\tdata2\n")

        output_file = os.path.join(temp_dir, "output.txt")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--column",
            "umi",
            "--sep",
            "\t",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            aligned_lines = [line.strip() for line in f.readlines() if line.strip()]
        assert len(aligned_lines) == 2


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
