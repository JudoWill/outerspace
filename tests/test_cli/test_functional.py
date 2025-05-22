"""Functional tests for the OUTERSPACE CLI workflow"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import shutil
from pathlib import Path
from argparse import Namespace
from outerspace.cli.commands.findseq import FindSeqCommand
from outerspace.cli.commands.collapse import CollapseCommand
from outerspace.cli.commands.count import CountCommand

@pytest.fixture
def temp_workspace():
    """Create a temporary workspace with test data structure"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create directory structure
        os.makedirs(os.path.join(temp_dir, 'reads/'))
        os.makedirs(os.path.join(temp_dir, 'results/extracted'))
        os.makedirs(os.path.join(temp_dir, 'results/collapsed'))
        os.makedirs(os.path.join(temp_dir, 'results/counted'))
        
        # Copy the real data files to temp directory
        # Note: These paths should be relative to the test file location
        data_files = [
            ('tests/data/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R1_001.fastq.gz'),
            ('tests/data/409-4_S1_L002_R2_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz'),
            ('tests/data/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz'),
            ('tests/data/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz'),
            ('tests/data/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz'),
            ('tests/data/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz'),
            ('tests/data/library_protospacers.txt', 'library_protospacers.txt')
        ]
        
        for src, dst in data_files:
            src_path = src
            dst_path = os.path.join(temp_dir, dst)
            if os.path.exists(src_path):
                shutil.copy2(src_path, dst_path)
            else:
                # Create empty files for testing if real files don't exist
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                open(dst_path, 'w').close()
        
        # Copy the config file
        config_src = 'tests/configs/grnaquery.cfg'
        config_dst = os.path.join(temp_dir, 'grnaquery.cfg')
        shutil.copy2(config_src, config_dst)

        yield temp_dir

def test_full_workflow(temp_workspace):
    """Test the full workflow from example.sh"""
    # Step 1: Run findseq commands for each sample
    
    pairs = [('reads/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz', 'shuffle'),
             ('reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'M1-lib'),
             ('reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'M2-lib')]
    
    for read1, read2, output_name in pairs:
        findseq_args = Namespace(
            command='findseq',
            config_filename=os.path.join(temp_workspace, 'grnaquery.cfg'),
            read1_filename=os.path.join(temp_workspace, read1),
            read2_filename=os.path.join(temp_workspace, read2),
            output_filename=os.path.join(temp_workspace, f'results/extracted/{output_name}.csv')
        )
        findseq_cmd = FindSeqCommand(findseq_args)
        findseq_cmd.run()
        
        # Verify findseq output
        assert os.path.exists(findseq_args.output_filename)
    
    # Step 2: Run collapse command
    
    collapse_args = Namespace(
                command='collapse',
                input_dir=os.path.join(temp_workspace, 'results/extracted'),
                output_dir=os.path.join(temp_workspace, 'results/collapsed'),
                columns='UMI_5prime,UMI_3prime',
                mismatches=2,
                row_limit=None,
                method='directional',
                sep=','
    )
    collapse_cmd = CollapseCommand(collapse_args)
    collapse_cmd.run()
    
    # Verify collapse output
    collapsed_files = os.listdir(collapse_args.output_dir)
    assert len(collapsed_files) > 0
    
    # Step 3: Run count command
    count_args = Namespace(
        command='count',
        input_dir=os.path.join(temp_workspace, 'results/collapsed'),
        output_dir=os.path.join(temp_workspace, 'results/counted'),
        barcode_column='UMI_5prime_UMI_3prime_corrected',
        key_column='protospacer',
        metrics=os.path.join(temp_workspace, 'results/counted/counts.yaml'),
        allowed_list=None,
        row_limit=None,
        downsample=None,
        detailed=False,
        sep=','
    )
    count_cmd = CountCommand(count_args)
    count_cmd.run()
    
    # Verify count output
    assert os.path.exists(count_args.metrics)
    counted_files = os.listdir(count_args.output_dir)
    assert len(counted_files) > 0

def test_workflow_with_allowed_list(temp_workspace):
    """Test the workflow with an allowed list for counting"""
    # Create allowed list file
    allowed_list_path = os.path.join(temp_workspace, 'library_protospacers.txt')
    
    
    # Run the workflow steps
    # Step 1: findseq (same as before)
    pairs = [('reads/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz', 'shuffle'),
             ('reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'M1-lib'),
             ('reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'M2-lib')]
    
    for read1, read2, output_name in pairs:
        findseq_args = Namespace(
            command='findseq',
            config_filename=os.path.join(temp_workspace, 'grnaquery.cfg'),
            read1_filename=os.path.join(temp_workspace, read1),
            read2_filename=os.path.join(temp_workspace, read2),
            output_filename=os.path.join(temp_workspace, f'results/extracted/{output_name}.csv')
        )
        findseq_cmd = FindSeqCommand(findseq_args)
        findseq_cmd.run()
    
    # Step 2: collapse (same as before)
    collapse_args = Namespace(
                command='collapse',
                input_dir=os.path.join(temp_workspace, 'results/extracted'),
                output_dir=os.path.join(temp_workspace, 'results/collapsed'),
                columns='UMI_5prime,UMI_3prime',
                mismatches=2,
                row_limit=None,
                method='directional',
                sep=','
    )
    collapse_cmd = CollapseCommand(collapse_args)
    collapse_cmd.run()
    
    # Step 3: count with allowed list
    count_args = Namespace(
        command='count',
        input_dir=os.path.join(temp_workspace, 'results/collapsed'),
        output_dir=os.path.join(temp_workspace, 'results/counted'),
        barcode_column='UMI_5prime_UMI_3prime_corrected',
        key_column='protospacer',
        metrics=os.path.join(temp_workspace, 'results/counted/counts.yaml'),
        allowed_list=allowed_list_path,
        row_limit=None,
        downsample=None,
        detailed=False,
        sep=','
    )
    
    count_cmd = CountCommand(count_args)
    count_cmd.run()
    
    # Verify count output with allowed list
    assert os.path.exists(count_args.metrics)
    counted_files = os.listdir(count_args.output_dir)
    assert len(counted_files) > 0

