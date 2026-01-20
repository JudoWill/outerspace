"""Tests for the align command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
import csv
from outerspace.cli.main import Cli


def test_align_initialization():
    """Test that align command initializes correctly"""
    args = [
        "align",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
    ]
    cli = Cli(args)
    assert cli.args.input_dir == "test_input"
    assert cli.args.output_dir == "test_output"
    assert cli.args.key_column == "sequence"
    assert cli.args.barcode_column == "barcode"
    assert cli.args.sep == ","
    assert cli.args.match == 5
    assert cli.args.mismatch == -4
    assert cli.args.gap == -8
    assert cli.args.algorithm == 1
    assert cli.args.min_count == 0
    assert cli.args.top_n is None
    assert cli.args.min_frequency == 0.0
    assert cli.args.align_by_barcode is False


def test_align_single_file_initialization():
    """Test that align command initializes correctly in single file mode"""
    args = [
        "align",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
    ]
    cli = Cli(args)
    assert cli.args.input_file == "test_input.csv"
    assert cli.args.output_file == "test_output.csv"
    assert cli.args.key_column == "sequence"
    assert cli.args.barcode_column == "barcode"
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
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_align_missing_key_column():
    """Test that align command handles missing key-column"""
    args = [
        "align",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--barcode-column",
        "barcode",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Please provide --key-column"):
        cli.run()


def test_align_missing_barcode_column():
    """Test that align command handles missing barcode-column"""
    args = [
        "align",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--key-column",
        "sequence",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Please provide --barcode-column"):
        cli.run()


def test_align_missing_output_dir_for_directory_mode():
    """Test that align command requires output-dir for directory mode"""
    args = [
        "align",
        "--input-dir",
        "test_input",
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
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
        "test_output.csv",
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
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
        "--key-column",
        "sequence",
        "--barcode-column",
        "barcode",
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
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            os.path.join(temp_dir, "output.csv"),
            "--key-column",
            "nonexistent_column",
            "--barcode-column",
            "header2",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Columns not found in input file"):
            cli.run()


def test_align_single_file_basic():
    """Test basic single file alignment functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with sequences and barcodes
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")
            f.write("AACTATA,BC3,data3\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)

        # Read and verify CSV output
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Should have 3 rows (one per unique sequence)
        assert len(rows) == 3

        # Verify CSV structure
        assert "sequence" in rows[0]
        assert "aligned_sequence" in rows[0]
        assert "unique_barcode_count" in rows[0]

        # All aligned sequences should have the same length
        aligned_lengths = set(len(row["aligned_sequence"]) for row in rows)
        assert len(aligned_lengths) == 1, "All aligned sequences should have the same length"


def test_align_single_file_stdout():
    """Test single file alignment output to stdout"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")

        args = [
            "align",
            "--input-file",
            input_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")
            f.write("AACTATA,BC3,data3\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 3


def test_align_directory_mode():
    """Test directory batch processing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create input directory with CSV files
        input_dir = os.path.join(temp_dir, "input")
        os.makedirs(input_dir)

        # Create first CSV file
        file1 = os.path.join(input_dir, "file1.csv")
        with open(file1, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")

        # Create second CSV file
        file2 = os.path.join(input_dir, "file2.csv")
        with open(file2, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("GGGTTTAA,BC3,data3\n")
            f.write("GGGTTTAG,BC4,data4\n")

        # Create output directory
        output_dir = os.path.join(temp_dir, "output")

        args = [
            "align",
            "--input-dir",
            input_dir,
            "--output-dir",
            output_dir,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            reader = csv.DictReader(f)
            rows1 = list(reader)
            assert len(rows1) == 2

        with open(output_file2, "r") as f:
            reader = csv.DictReader(f)
            rows2 = list(reader)
            assert len(rows2) == 2


def test_align_empty_sequences():
    """Test alignment with empty sequences in CSV"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with some empty sequences
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write(",BC2,data2\n")  # Empty sequence
            f.write("AACTATA,BC3,data3\n")
            f.write(",BC4,data4\n")  # Empty sequence

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
        ]
        cli = Cli(args)
        cli.run()

        # Should only align non-empty sequences
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        # Should have 2 rows (empty ones filtered out)
        assert len(rows) == 2


def test_align_single_sequence():
    """Test alignment with only one sequence"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with single sequence
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
        ]
        cli = Cli(args)
        cli.run()

        # Single sequence should be returned as-is
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["sequence"] == "AACTTATA"
        assert rows[0]["unique_barcode_count"] == "1"


def test_align_row_limit():
    """Test alignment with row limit"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")
            f.write("AACTATA,BC3,data3\n")
            f.write("GGGTTTAA,BC4,data4\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--row-limit",
            "3",  # Only process first 3 data rows
        ]
        cli = Cli(args)
        cli.run()

        # Should only align sequences from first 3 data rows
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 3


def test_align_different_algorithms():
    """Test alignment with different algorithm types"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATG,BC2,data2\n")
            f.write("AACTATA,BC3,data3\n")

        # Test local algorithm (0)
        output_file1 = os.path.join(temp_dir, "output_local.csv")
        args1 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file1,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--algorithm",
            "0",
        ]
        cli1 = Cli(args1)
        cli1.run()

        # Test global algorithm (1)
        output_file2 = os.path.join(temp_dir, "output_global.csv")
        args2 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file2,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--algorithm",
            "1",
        ]
        cli2 = Cli(args2)
        cli2.run()

        # Test semi-global algorithm (2)
        output_file3 = os.path.join(temp_dir, "output_semiglobal.csv")
        args3 = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file3,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-dir",
            os.path.join(temp_dir, "output"),
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
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
            f.write("sequence\tbarcode\tother\n")
            f.write("AACTTATA\tBC1\tdata1\n")
            f.write("AACTTATG\tBC2\tdata2\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--sep",
            "\t",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2


def test_align_min_count_filter():
    """Test alignment with min-count filtering"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with sequences having different barcode counts
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            # Sequence1 appears with 3 barcodes
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATA,BC2,data2\n")
            f.write("AACTTATA,BC3,data3\n")
            # Sequence2 appears with 1 barcode
            f.write("GGGTTTAA,BC4,data4\n")
            # Sequence3 appears with 2 barcodes
            f.write("CCCTTTGG,BC5,data5\n")
            f.write("CCCTTTGG,BC6,data6\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--min-count",
            "2",  # Only sequences with at least 2 unique barcodes
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        # Should have 2 sequences (AACTTATA with 3 barcodes, CCCTTTGG with 2 barcodes)
        assert len(rows) == 2
        sequences = {row["sequence"] for row in rows}
        assert "AACTTATA" in sequences
        assert "CCCTTTGG" in sequences
        assert "GGGTTTAA" not in sequences


def test_align_top_n_filter():
    """Test alignment with top-n filtering"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with sequences having different barcode counts
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            # Sequence1 appears with 3 barcodes
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATA,BC2,data2\n")
            f.write("AACTTATA,BC3,data3\n")
            # Sequence2 appears with 1 barcode
            f.write("GGGTTTAA,BC4,data4\n")
            # Sequence3 appears with 2 barcodes
            f.write("CCCTTTGG,BC5,data5\n")
            f.write("CCCTTTGG,BC6,data6\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--top-n",
            "2",  # Top 2 sequences by barcode count
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        # Should have 2 sequences (AACTTATA with 3, CCCTTTGG with 2)
        assert len(rows) == 2
        sequences = {row["sequence"] for row in rows}
        assert "AACTTATA" in sequences
        assert "CCCTTTGG" in sequences
        assert "GGGTTTAA" not in sequences


def test_align_min_frequency_filter():
    """Test alignment with min-frequency filtering"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            # Sequence1 appears with 3 barcodes (50% of total 6 unique barcodes)
            f.write("AACTTATA,BC1,data1\n")
            f.write("AACTTATA,BC2,data2\n")
            f.write("AACTTATA,BC3,data3\n")
            # Sequence2 appears with 1 barcode (16.7%)
            f.write("GGGTTTAA,BC4,data4\n")
            # Sequence3 appears with 2 barcodes (33.3%)
            f.write("CCCTTTGG,BC5,data5\n")
            f.write("CCCTTTGG,BC6,data6\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--min-frequency",
            "30.0",  # At least 30% frequency
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        # Should have 2 sequences (AACTTATA 50%, CCCTTTGG 33.3%)
        assert len(rows) == 2
        sequences = {row["sequence"] for row in rows}
        assert "AACTTATA" in sequences
        assert "CCCTTTGG" in sequences
        assert "GGGTTTAA" not in sequences


def test_align_by_barcode():
    """Test alignment with --align-by-barcode option"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            # Sequence1 appears with barcode BC1
            f.write("AACTTATA,BC1,data1\n")
            # Sequence2 appears with barcode BC1 (same barcode as sequence1)
            f.write("GGGTTTAA,BC1,data2\n")
            # Sequence1 also appears with barcode BC2 (multiple barcodes)
            f.write("AACTTATA,BC2,data3\n")
            # Sequence3 appears with barcode BC3
            f.write("CCCTTTGG,BC3,data4\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--align-by-barcode",
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        # Should have 4 rows:
        # - AACTTATA with BC1
        # - GGGTTTAA with BC1
        # - AACTTATA with BC2 (appears again because it has multiple barcodes)
        # - CCCTTTGG with BC3
        assert len(rows) == 4
        
        # Verify CSV structure includes barcode column
        assert "barcode" in rows[0]
        assert "sequence" in rows[0]
        assert "aligned_sequence" in rows[0]
        assert "unique_barcode_count" in rows[0]
        
        # Verify AACTTATA appears twice (once per barcode)
        aacttata_rows = [row for row in rows if row["sequence"] == "AACTTATA"]
        assert len(aacttata_rows) == 2
        barcodes = {row["barcode"] for row in aacttata_rows}
        assert barcodes == {"BC1", "BC2"}


def test_align_multiple_filters():
    """Test alignment with multiple filters applied simultaneously"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            # Sequence1 appears with 5 barcodes
            for i in range(5):
                f.write(f"AACTTATA,BC{i+1},data{i+1}\n")
            # Sequence2 appears with 3 barcodes
            for i in range(3):
                f.write(f"GGGTTTAA,BC{i+6},data{i+6}\n")
            # Sequence3 appears with 2 barcodes
            for i in range(2):
                f.write(f"CCCTTTGG,BC{i+9},data{i+9}\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--min-count",
            "2",  # At least 2 barcodes
            "--top-n",
            "2",  # Top 2 sequences
            "--min-frequency",
            "20.0",  # At least 20% frequency
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        # Should have 2 sequences (AACTTATA with 5, GGGTTTAA with 3)
        # Both pass min-count (>=2), both are in top-2, both pass min-frequency (>=20%)
        assert len(rows) == 2
        sequences = {row["sequence"] for row in rows}
        assert "AACTTATA" in sequences
        assert "GGGTTTAA" in sequences
        assert "CCCTTTGG" not in sequences  # Filtered out by top-n


def test_align_empty_after_filtering():
    """Test alignment when all sequences are filtered out"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("sequence,barcode,other\n")
            f.write("AACTTATA,BC1,data1\n")
            f.write("GGGTTTAA,BC2,data2\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "align",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--key-column",
            "sequence",
            "--barcode-column",
            "barcode",
            "--min-count",
            "5",  # Both sequences have only 1 barcode, so both filtered out
        ]
        cli = Cli(args)
        cli.run()

        # Should create empty output file (just header)
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 0


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
