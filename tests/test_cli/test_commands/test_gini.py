"""Tests for the gini command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from argparse import Namespace
from outerspace.cli.commands.gini import GiniCommand

def test_gini_initialization():
    """Test that gini command initializes correctly"""
    args = Namespace(
        command='gini',
        input_file='test_input.csv',
        column='counts',
        count_column=None,
        scale=None,
        sep=',',
        allowed_list=None
    )
    cmd = GiniCommand(args)
    assert cmd.args == args

def test_gini_missing_input_file():
    """Test that gini command handles missing input file"""
    args = Namespace(
        command='gini',
        input_file=None,
        column='counts',
        count_column=None,
        scale=None,
        sep=',',
        allowed_list=None
    )
    cmd = GiniCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_gini_missing_column():
    """Test that gini command handles missing column"""
    args = Namespace(
        command='gini',
        input_file='test_input.csv',
        column=None,
        count_column=None,
        scale=None,
        sep=',',
        allowed_list=None
    )
    cmd = GiniCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_gini_invalid_scale():
    """Test that gini command handles invalid scale value"""
    args = Namespace(
        command='gini',
        input_file='test_input.csv',
        column='counts',
        count_column=None,
        scale=-1,
        sep=',',
        allowed_list=None
    )
    cmd = GiniCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_gini_nonexistent_input_file():
    """Test that gini command handles nonexistent input file"""
    args = Namespace(
        command='gini',
        input_file='nonexistent.csv',
        column='counts',
        count_column=None,
        scale=None,
        sep=',',
        allowed_list=None
    )
    cmd = GiniCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_gini_nonexistent_column():
    """Test that gini command handles nonexistent column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = Namespace(
            command='gini',
            input_file=temp_file.name,
            column='nonexistent',
            count_column=None,
            scale=None,
            sep=',',
            allowed_list=None
        )
        cmd = GiniCommand(args)
        with pytest.raises(ValueError):
            cmd.run()

def test_gini_nonexistent_count_column():
    """Test that gini command handles nonexistent count column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = Namespace(
            command='gini',
            input_file=temp_file.name,
            column='header1',
            count_column='nonexistent',
            scale=None,
            sep=',',
            allowed_list=None
        )
        cmd = GiniCommand(args)
        with pytest.raises(ValueError):
            cmd.run()

def test_gini_nonexistent_allowed_list():
    """Test that gini command handles nonexistent allowed list file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = Namespace(
            command='gini',
            input_file=temp_file.name,
            column='header1',
            count_column=None,
            scale=None,
            sep=',',
            allowed_list='nonexistent.txt'
        )
        cmd = GiniCommand(args)
        with pytest.raises(ValueError):
            cmd.run()
