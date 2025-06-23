"""Tests for the visualize command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.cli.main import Cli


@pytest.mark.skip(reason="Not implemented")
def test_visualize_initialization():
    """Test that visualize command initializes correctly"""
    args = ["visualize", "test_input", "test_output", "--bins", "50", "--format", "png"]
    cli = Cli(args)
    assert cli.args.input_dir == "test_input"
    assert cli.args.output_dir == "test_output"
    assert cli.args.bins == 50
    assert cli.args.format == "png"
    assert cli.args.sep == ","
    assert cli.args.title_prefix is None
    assert cli.args.xlabel is None
    assert cli.args.ylabel is None
    assert cli.args.log_scale is False


@pytest.mark.skip(reason="Not implemented")
def test_visualize_missing_input_dir():
    """Test that visualize command handles missing input directory"""
    args = [
        "visualize",
        "--output-dir",
        "test_output",
        "--bins",
        "50",
        "--format",
        "png",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


@pytest.mark.skip(reason="Not implemented")
def test_visualize_missing_output_dir():
    """Test that visualize command handles missing output directory"""
    args = ["visualize", "test_input", "--bins", "50", "--format", "png"]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


@pytest.mark.skip(reason="Not implemented")
def test_visualize_invalid_bins():
    """Test that visualize command handles invalid number of bins"""
    args = ["visualize", "test_input", "test_output", "--bins", "0", "--format", "png"]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


@pytest.mark.skip(reason="Not implemented")
def test_visualize_invalid_format():
    """Test that visualize command handles invalid output format"""
    args = [
        "visualize",
        "test_input",
        "test_output",
        "--bins",
        "50",
        "--format",
        "invalid_format",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_visualize_not_implemented():
    """Test that visualize command is not yet implemented"""
    args = ["visualize", "test_input", "test_output", "--bins", "50", "--format", "png"]
    cli = Cli(args)
    with pytest.raises(NotImplementedError):
        cli.run()
