"""Tests for the gini command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli

def test_gini_initialization():
    """Test that gini command initializes correctly"""
    args = [
        'gini',
        'test_input.csv',
        '--column', 'counts'
    ]
    cli = Cli(args)
    assert cli.args.input_file == 'test_input.csv'
    assert cli.args.column == 'counts'
    assert cli.args.sep == ','
    assert cli.args.scale is None
    assert cli.args.count_column is None
    assert cli.args.allowed_list is None

def test_gini_missing_input_file():
    """Test that gini command handles missing input file"""
    args = [
        'gini',
        '--column', 'counts'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_gini_missing_column():
    """Test that gini command handles missing column"""
    args = [
        'gini',
        'test_input.csv'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_gini_invalid_scale():
    """Test that gini command handles invalid scale value"""
    args = [
        'gini',
        'test_input.csv',
        '--column', 'counts',
        '--scale', '-1'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_gini_nonexistent_input_file():
    """Test that gini command handles nonexistent input file"""
    args = [
        'gini',
        'nonexistent.csv',
        '--column', 'counts'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_gini_nonexistent_column():
    """Test that gini command handles nonexistent column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_gini_nonexistent_count_column():
    """Test that gini command handles nonexistent count column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'header1',
            '--count-column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_gini_nonexistent_allowed_list():
    """Test that gini command handles nonexistent allowed list file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'header1',
            '--allowed-list', 'nonexistent.txt'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()
