"""Tests for the collapse command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.commands.collapse import CollapseCommand

@pytest.mark.skip(reason="Not implemented")
def test_collapse_initialization():
    """Test that collapse command initializes correctly"""
    args = Namespace(
        command='collapse',
        input_dir='test_input',
        output_dir='test_output',
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    assert cmd.args == args

@pytest.mark.skip(reason="Not implemented")
def test_collapse_missing_input_dir(capsys):
    """Test that collapse command handles missing input directory"""
    args = Namespace(
        command='collapse',
        input_dir=None,
        output_dir='test_output',
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide an input directory" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_collapse_missing_output_dir(capsys):
    """Test that collapse command handles missing output directory"""
    args = Namespace(
        command='collapse',
        input_dir='test_input',
        output_dir=None,
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide an output directory" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_collapse_missing_columns(capsys):
    """Test that collapse command handles missing columns"""
    args = Namespace(
        command='collapse',
        input_dir='test_input',
        output_dir='test_output',
        columns=None,
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide columns to correct" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_collapse_invalid_method(capsys):
    """Test that collapse command handles invalid clustering method"""
    args = Namespace(
        command='collapse',
        input_dir='test_input',
        output_dir='test_output',
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='invalid_method'
    )
    cmd = CollapseCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Invalid clustering method" in captured.out

def test_collapse_not_implemented():
    """Test that collapse command is not yet implemented"""
    args = Namespace(
        command='collapse',
        input_dir='test_input',
        output_dir='test_output',
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    with pytest.raises(NotImplementedError):
        cmd.run() 