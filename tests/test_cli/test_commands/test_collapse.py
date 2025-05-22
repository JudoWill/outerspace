"""Tests for the collapse command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from argparse import Namespace
from outerspace.cli.commands.collapse import CollapseCommand

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

def test_collapse_missing_input_dir():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_collapse_missing_output_dir():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_collapse_missing_columns():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_collapse_invalid_method():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_collapse_nonexistent_input_dir():
    """Test that collapse command handles nonexistent input directory"""
    args = Namespace(
        command='collapse',
        input_dir='nonexistent_dir',
        output_dir='test_output',
        columns='umi3,umi5',
        mismatches=2,
        sep=',',
        row_limit=None,
        method='directional'
    )
    cmd = CollapseCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_collapse_empty_input_dir():
    """Test that collapse command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = Namespace(
            command='collapse',
            input_dir=temp_dir,
            output_dir='test_output',
            columns='umi3,umi5',
            mismatches=2,
            sep=',',
            row_limit=None,
            method='directional'
        )
        cmd = CollapseCommand(args)
        with pytest.raises(ValueError):
            cmd.run()

def test_collapse_parse_columns():
    """Test that column parsing works correctly"""
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
    columns = cmd._parse_columns(args.columns)
    assert columns == ['umi3', 'umi5']
    
    # Test with spaces
    columns = cmd._parse_columns('umi3, umi5')
    assert columns == ['umi3', 'umi5']
    
    # Test with single column
    columns = cmd._parse_columns('umi3')
    assert columns == ['umi3']
    