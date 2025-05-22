"""Tests for the count command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.commands.count import CountCommand

def test_count_initialization():
    """Test that count command initializes correctly"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir='test_output',
        barcode_column='umi3_umi5_corrected',
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    assert cmd.args == args

@pytest.mark.skip(reason="Not implemented")
def test_count_missing_input_dir(capsys):
    """Test that count command handles missing input directory"""
    args = Namespace(
        command='count',
        input_dir=None,
        output_dir='test_output',
        barcode_column='umi3_umi5_corrected',
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide an input directory" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_count_missing_output_dir(capsys):
    """Test that count command handles missing output directory"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir=None,
        barcode_column='umi3_umi5_corrected',
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide an output directory" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_count_missing_barcode_column(capsys):
    """Test that count command handles missing barcode column"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir='test_output',
        barcode_column=None,
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide a barcode column" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_count_missing_key_column(capsys):
    """Test that count command handles missing key column"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir='test_output',
        barcode_column='umi3_umi5_corrected',
        key_column=None,
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Please provide a key column" in captured.out

@pytest.mark.skip(reason="Not implemented")
def test_count_invalid_downsample(capsys):
    """Test that count command handles invalid downsample value"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir='test_output',
        barcode_column='umi3_umi5_corrected',
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=-1,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    cmd.run()
    captured = capsys.readouterr()
    assert "Downsample value must be positive" in captured.out

def test_count_not_implemented():
    """Test that count command is not yet implemented"""
    args = Namespace(
        command='count',
        input_dir='test_input',
        output_dir='test_output',
        barcode_column='umi3_umi5_corrected',
        key_column='protospacer',
        sep=',',
        row_limit=None,
        allowed_list=None,
        detailed=False,
        downsample=None,
        random_seed=None,
        metrics=None
    )
    cmd = CountCommand(args)
    with pytest.raises(NotImplementedError):
        cmd.run() 