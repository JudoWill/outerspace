"""Tests for the stats command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import csv
from outerspace.cli.main import Cli


def test_stats_initialization():
    """Test that stats command initializes correctly with both config and CLI args"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file based on grnaquery.toml
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")
            f.write("B,2\n")
            f.write("C,3\n")

        args = ["stats", "--config", config_file, input_file]
        cli = Cli(args)
        cli.run()


def test_stats_cli_override_config():
    """Test that CLI arguments override config file settings"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        # Create test input file with different column names
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            f.write("A,1\n")
            f.write("B,2\n")
            f.write("C,3\n")

        args = ["stats", input_file, "--key-column", "key", "--count-column", "count"]
        cli = Cli(args)
        cli.run()


def test_stats_missing_required_args():
    """Test that stats command requires either config or CLI args"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            f.write("A,1\n")

        args = ["stats", input_file]
        with pytest.raises(
            ValueError, match="Please provide either --key-column or --config"
        ):
            cli = Cli(args)
            cli.run()


def test_stats_invalid_scale():
    """Test that stats command handles invalid scale value"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")

        args = ["stats", "--config", config_file, "--scale", "-1", input_file]
        with pytest.raises(ValueError, match="Scale value must be positive"):
            cli = Cli(args)
            cli.run()


def test_stats_nonexistent_input_file():
    """Test that stats command handles nonexistent input file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        args = ["stats", "--config", config_file, "nonexistent.csv"]
        with pytest.raises(ValueError, match="No files were successfully processed"):
            cli = Cli(args)
            cli.run()


def test_stats_nonexistent_columns():
    """Test that stats command handles nonexistent columns"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")

        args = [
            "stats",
            input_file,
            "--key-column",
            "nonexistent",
            "--count-column",
            "nonexistent",
        ]
        with pytest.raises(ValueError):
            cli = Cli(args)
            cli.run()


def test_stats_with_allowed_list():
    """Test that stats command works with allowed list"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("A\nB\n")

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")
            f.write("B,2\n")
            f.write("C,3\n")

        args = [
            "stats",
            "--config",
            config_file,
            "--allowed-list",
            allowed_list,
            input_file,
        ]
        cli = Cli(args)
        cli.run()


def test_stats_multiple_input_files():
    """Test that stats command handles multiple input files"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[stats]
key_column = "protospacer"
count_column = "UMI_5prime_UMI_3prime_corrected_count"
"""
            )

        # Create test input files
        input_files = []
        for i in range(3):
            input_file = os.path.join(temp_dir, f"test{i}.csv")
            with open(input_file, "w") as f:
                f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
                f.write("A,1\n")
                f.write("B,2\n")
                f.write("C,3\n")
            input_files.append(input_file)

        args = ["stats", "--config", config_file] + input_files
        cli = Cli(args)
        cli.run()
