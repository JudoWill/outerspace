"""Tests for the findseq command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
import csv
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

def test_findseq_with_example_data():
    """Test findseq command with real example data"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        args = Namespace(
            command='findseq',
            config_filename='tests/configs/grnaquery.cfg',
            read1_filename='tests/data/409-4_S1_L002_R1_001.fastq.gz',
            read2_filename='tests/data/409-4_S1_L002_R2_001.fastq.gz',
            output_filename=os.path.join(temp_dir, 'shuffle.csv'),
            fastqfiles=None,
            outdir=None
        )
        cmd = FindSeqCommand(args)
        cmd.run()
    
        # Verify output file exists and has content
        assert os.path.exists(os.path.join(temp_dir, 'shuffle.csv'))
        
        with open(os.path.join(temp_dir, 'shuffle.csv'), 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == ['read_id', 'UMI_5prime', 'protospacer', 'downstreamof_protospacer', 'UMI_3prime']
            assert len(list(reader)) == 442
        
        
