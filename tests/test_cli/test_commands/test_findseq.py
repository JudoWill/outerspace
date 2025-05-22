"""Tests for the findseq command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from argparse import Namespace
from outerspace.cli.commands.findseq import FindSeqCommand

def test_findseq_initialization():
    """Test that findseq command initializes correctly"""
    args = Namespace(
        command='findseq',
        config_filename='test_config.toml',
        read1_filename='test_r1.fastq',
        read2_filename='test_r2.fastq',
        output_filename='test_output.csv',
        fastqfiles=None,
        outdir=None
    )
    cmd = FindSeqCommand(args)
    assert cmd.args == args

def test_findseq_missing_config(capsys):
    """Test that findseq command handles missing config file"""
    args = Namespace(
        command='findseq',
        config_filename=None,
        read1_filename='test_r1.fastq',
        read2_filename='test_r2.fastq',
        output_filename='test_output.csv',
        fastqfiles=None,
        outdir=None
    )
    cmd = FindSeqCommand(args)
    with pytest.raises(ValueError):
        cmd.run()

def test_findseq_missing_reads(capsys):
    """Test that findseq command handles missing read files"""
    args = Namespace(
        command='findseq',
        config_filename='test_config.toml',
        read1_filename=None,
        read2_filename=None,
        output_filename=None,
        fastqfiles=None,
        outdir=None
    )
    cmd = FindSeqCommand(args)  
    with pytest.raises(ValueError):
        cmd.run()

def test_findseq_single_read_not_implemented():
    """Test that single read processing is not implemented"""
    args = Namespace(
        command='findseq',
        config_filename='test_config.toml',
        read1_filename='test_r1.fastq',
        read2_filename=None,
        output_filename='test_output.csv',
        fastqfiles=None,
        outdir=None
    )
    cmd = FindSeqCommand(args)
    with pytest.raises(NotImplementedError):
        cmd.run()

def test_findseq_no_reads(capsys):
    """Test that findseq command handles no reads case"""
    args = Namespace(
        command='findseq',
        config_filename='test_config.toml',
        read1_filename=None,
        read2_filename=None,
        output_filename='test_output.csv',
        fastqfiles=None,
        outdir=None
    )
    cmd = FindSeqCommand(args)
    with pytest.raises(ValueError):
        cmd.run()
