"""Tests for the count command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
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

def test_count_missing_input_dir():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_missing_output_dir():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_missing_barcode_column():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_missing_key_column():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_invalid_downsample():
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_nonexistent_input_dir():
    """Test that count command handles nonexistent input directory"""
    args = Namespace(
        command='count',
        input_dir='nonexistent_dir',
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
    with pytest.raises(ValueError):
        cmd.run()

def test_count_empty_input_dir():
    """Test that count command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = Namespace(
            command='count',
            input_dir=temp_dir,
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
        with pytest.raises(ValueError):
            cmd.run()

