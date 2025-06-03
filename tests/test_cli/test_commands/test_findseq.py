"""Tests for the findseq command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
import csv
from outerspace.cli.main import Cli 

def test_findseq_initialization():
    """Test that findseq command initializes correctly"""
    args = [
        'findseq',
        'test_config.toml',
        '-1', 'test_r1.fastq',
        '-2', 'test_r2.fastq',
        '-o', 'test_output.csv'
    ]
    cli = Cli(args)
    assert cli.args.config == 'test_config.toml'
    assert cli.args.read1_filename == 'test_r1.fastq'
    assert cli.args.read2_filename == 'test_r2.fastq'
    assert cli.args.output_filename == 'test_output.csv'

def test_findseq_missing_reads():
    """Test that findseq command handles missing read files"""
    args = [
        'findseq',
        'test_config.toml',
        '-o', 'test_output.csv'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_findseq_no_reads():
    """Test that findseq command handles no reads case"""
    args = [
        'findseq',
        'test_config.toml',
        '-o', 'test_output.csv'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_findseq_with_example_data():
    """Test findseq command with real example data"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            'findseq',
            'tests/configs/grnaquery.cfg',
            '-1', 'tests/data/409-4_S1_L002_R1_001.fastq.gz',
            '-2', 'tests/data/409-4_S1_L002_R2_001.fastq.gz',
            '-o', os.path.join(temp_dir, 'shuffle.csv')
        ]
        cli = Cli(args)
        cli.run()
    
        # Verify output file exists and has content
        assert os.path.exists(os.path.join(temp_dir, 'shuffle.csv'))
        
        with open(os.path.join(temp_dir, 'shuffle.csv'), 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == ['read_id', 'UMI_5prime', 'protospacer', 'downstreamof_protospacer', 'UMI_3prime']
            assert len(list(reader)) == 442
        
        
