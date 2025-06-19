"""Functional tests for the OUTERSPACE CLI workflow using iterative commands"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
from outerspace.cli.main import Cli

def test_full_workflow(temp_workspace):
    """Test the full workflow from example.sh"""
    # Step 1: Run findseq commands for each sample
    pairs = [('reads/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz', 'shuffle'),
             ('reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'M1-lib'),
             ('reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'M2-lib')]

    for read1, read2, output_name in pairs:
        findseq_args = [
            'findseq',
            os.path.join(temp_workspace, 'grnaquery.toml'),
            '-1', os.path.join(temp_workspace, read1),
            '-2', os.path.join(temp_workspace, read2),
            '-o', os.path.join(temp_workspace, f'findseq/{output_name}.csv')
        ]
        cli = Cli(findseq_args)
        cli.run()

        # Verify findseq output
        assert os.path.exists(os.path.join(temp_workspace, f'findseq/{output_name}.csv'))

    # Step 2: Run collapse command
    collapse_args = [
        'collapse',
        '--input-dir', os.path.join(temp_workspace, 'findseq'),
        '--output-dir', os.path.join(temp_workspace, 'collapse'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'collapse/{output_name}.csv'))

    # Step 3: Run count command
    count_args = [
        'count',
        '--input-dir', os.path.join(temp_workspace, 'collapse'),
        '--output-dir', os.path.join(temp_workspace, 'count'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'count/{output_name}.csv'))

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
        findseq_args = [
            'findseq',
            os.path.join(temp_workspace, 'grnaquery.toml'),
            '-1', os.path.join(temp_workspace, read1),
            '-2', os.path.join(temp_workspace, read2),
            '-o', os.path.join(temp_workspace, f'findseq/{output_name}.csv')
        ]
        cli = Cli(findseq_args)
        cli.run()

    # Step 2: collapse (same as before)
    collapse_args = [
        'collapse',
        '--input-dir', os.path.join(temp_workspace, 'findseq'),
        '--output-dir', os.path.join(temp_workspace, 'collapse'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: count with allowed list
    count_args = [
        'count',
        '--input-dir', os.path.join(temp_workspace, 'collapse'),
        '--output-dir', os.path.join(temp_workspace, 'count'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--allowed-list', allowed_list_path
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'count/{output_name}.csv'))

def test_single_file_workflow(temp_workspace):
    """Test the workflow using single file mode for collapse and count"""
    # Step 1: Run findseq for a single sample
    findseq_args = [
        'findseq',
        os.path.join(temp_workspace, 'grnaquery.toml'),
        '-1', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R1_001.fastq.gz'),
        '-2', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R2_001.fastq.gz'),
        '-o', os.path.join(temp_workspace, 'findseq/shuffle.csv')
    ]
    cli = Cli(findseq_args)
    cli.run()

    # Verify findseq output
    assert os.path.exists(os.path.join(temp_workspace, 'findseq/shuffle.csv'))

    # Step 2: Run collapse command on single file
    collapse_args = [
        'collapse',
        '--input-file', os.path.join(temp_workspace, 'findseq/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'collapse/shuffle.csv'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional',
        '--metrics', os.path.join(temp_workspace, 'collapse/shuffle_metrics.yaml')
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    assert os.path.exists(os.path.join(temp_workspace, 'collapse/shuffle.csv'))
    assert os.path.exists(os.path.join(temp_workspace, 'collapse/shuffle_metrics.yaml'))

    # Step 3: Run count command on single file
    count_args = [
        'count',
        '--input-file', os.path.join(temp_workspace, 'collapse/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'count/shuffle.csv'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--metrics', os.path.join(temp_workspace, 'count/shuffle_metrics.yaml')
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    assert os.path.exists(os.path.join(temp_workspace, 'count/shuffle.csv'))

def test_single_file_workflow_with_allowed_list(temp_workspace):
    """Test the single file workflow with an allowed list for counting"""
    # Create allowed list file
    allowed_list_path = os.path.join(temp_workspace, 'library_protospacers.txt')

    # Step 1: Run findseq for a single sample
    findseq_args = [
        'findseq',
        os.path.join(temp_workspace, 'grnaquery.toml'),
        '-1', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R1_001.fastq.gz'),
        '-2', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R2_001.fastq.gz'),
        '-o', os.path.join(temp_workspace, 'findseq/shuffle.csv')
    ]
    cli = Cli(findseq_args)
    cli.run()

    # Step 2: Run collapse command on single file
    collapse_args = [
        'collapse',
        '--input-file', os.path.join(temp_workspace, 'findseq/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'collapse/shuffle.csv'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional',
        '--metrics', os.path.join(temp_workspace, 'collapse/shuffle_metrics.yaml')
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: Run count command on single file with allowed list
    count_args = [
        'count',
        '--input-file', os.path.join(temp_workspace, 'collapse/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'count/shuffle.csv'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--allowed-list', allowed_list_path,
        '--metrics', os.path.join(temp_workspace, 'count/shuffle_metrics.yaml')
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify outputs
    assert os.path.exists(os.path.join(temp_workspace, 'count/shuffle.csv'))
    assert os.path.exists(os.path.join(temp_workspace, 'count/shuffle_metrics.yaml')) 