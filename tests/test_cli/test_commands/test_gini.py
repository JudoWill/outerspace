"""Tests for the gini command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.commands.gini import GiniCommand

@pytest.mark.skip(reason="Not implemented")
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

@pytest.mark.skip(reason="Not implemented")
def test_gini_missing_input_file(capsys):
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
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide an input file" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_gini_missing_column(capsys):
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
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide a column to calculate Gini coefficient for" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_gini_invalid_scale(capsys):
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
    cmd.run()
    captured = capsys.readouterr()
    assert "Scale value must be positive" in captured.out

def test_gini_not_implemented():
    """Test that gini command is not yet implemented"""
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
    with pytest.raises(NotImplementedError):
        cmd.run() 