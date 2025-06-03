"""Functional tests for the OUTERSPACE CLI workflow"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import shutil
from pathlib import Path
from outerspace.cli.main import Cli

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
        findseq_args = [
            'findseq',
            os.path.join(temp_workspace, 'grnaquery.cfg'),
            '-1', os.path.join(temp_workspace, read1),
            '-2', os.path.join(temp_workspace, read2),
            '-o', os.path.join(temp_workspace, f'results/extracted/{output_name}.csv')
        ]
        cli = Cli(findseq_args)
        cli.run()

        # Verify findseq output
        assert os.path.exists(os.path.join(temp_workspace, f'results/extracted/{output_name}.csv'))

    # Step 2: Run collapse command
    collapse_args = [
        'collapse',
        '--input-dir', os.path.join(temp_workspace, 'results/extracted'),
        '--output-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'results/collapsed/{output_name}.csv'))

    # Step 3: Run count command
    count_args = [
        'count',
        '--input-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--output-dir', os.path.join(temp_workspace, 'results/counted'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'results/counted/{output_name}.csv'))

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
            os.path.join(temp_workspace, 'grnaquery.cfg'),
            '-1', os.path.join(temp_workspace, read1),
            '-2', os.path.join(temp_workspace, read2),
            '-o', os.path.join(temp_workspace, f'results/extracted/{output_name}.csv')
        ]
        cli = Cli(findseq_args)
        cli.run()

    # Step 2: collapse (same as before)
    collapse_args = [
        'collapse',
        '--input-dir', os.path.join(temp_workspace, 'results/extracted'),
        '--output-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: count with allowed list
    count_args = [
        'count',
        '--input-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--output-dir', os.path.join(temp_workspace, 'results/counted'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--allowed-list', allowed_list_path
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'results/counted/{output_name}.csv'))

def test_single_file_workflow(temp_workspace):
    """Test the workflow using single file mode for collapse and count"""
    # Step 1: Run findseq for a single sample
    findseq_args = [
        'findseq',
        os.path.join(temp_workspace, 'grnaquery.cfg'),
        '-1', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R1_001.fastq.gz'),
        '-2', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R2_001.fastq.gz'),
        '-o', os.path.join(temp_workspace, 'results/extracted/shuffle.csv')
    ]
    cli = Cli(findseq_args)
    cli.run()

    # Verify findseq output
    assert os.path.exists(os.path.join(temp_workspace, 'results/extracted/shuffle.csv'))

    # Step 2: Run collapse command on single file
    collapse_args = [
        'collapse',
        '--input-file', os.path.join(temp_workspace, 'results/extracted/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'results/collapsed/shuffle.csv'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional',
        '--metrics', os.path.join(temp_workspace, 'results/collapsed/shuffle_metrics.yaml')
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    assert os.path.exists(os.path.join(temp_workspace, 'results/collapsed/shuffle.csv'))
    assert os.path.exists(os.path.join(temp_workspace, 'results/collapsed/shuffle_metrics.yaml'))

    # Step 3: Run count command on single file
    count_args = [
        'count',
        '--input-file', os.path.join(temp_workspace, 'results/collapsed/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'results/counted/shuffle.csv'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--metrics', os.path.join(temp_workspace, 'results/counted/shuffle_metrics.yaml')
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    assert os.path.exists(os.path.join(temp_workspace, 'results/counted/shuffle.csv'))

def test_single_file_workflow_with_allowed_list(temp_workspace):
    """Test the single file workflow with an allowed list for counting"""
    # Create allowed list file
    allowed_list_path = os.path.join(temp_workspace, 'library_protospacers.txt')

    # Step 1: Run findseq for a single sample
    findseq_args = [
        'findseq',
        os.path.join(temp_workspace, 'grnaquery.cfg'),
        '-1', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R1_001.fastq.gz'),
        '-2', os.path.join(temp_workspace, 'reads/409-4_S1_L002_R2_001.fastq.gz'),
        '-o', os.path.join(temp_workspace, 'results/extracted/shuffle.csv')
    ]
    cli = Cli(findseq_args)
    cli.run()

    # Step 2: Run collapse command on single file
    collapse_args = [
        'collapse',
        '--input-file', os.path.join(temp_workspace, 'results/extracted/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'results/collapsed/shuffle.csv'),
        '--columns', 'UMI_5prime,UMI_3prime',
        '--mismatches', '2',
        '--method', 'directional',
        '--metrics', os.path.join(temp_workspace, 'results/collapsed/shuffle_metrics.yaml')
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: Run count command on single file with allowed list
    count_args = [
        'count',
        '--input-file', os.path.join(temp_workspace, 'results/collapsed/shuffle.csv'),
        '--output-file', os.path.join(temp_workspace, 'results/counted/shuffle.csv'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',
        '--key-column', 'protospacer',
        '--allowed-list', allowed_list_path,
        '--metrics', os.path.join(temp_workspace, 'results/counted/shuffle_metrics.yaml')
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify outputs
    assert os.path.exists(os.path.join(temp_workspace, 'results/counted/shuffle.csv'))
    assert os.path.exists(os.path.join(temp_workspace, 'results/counted/shuffle_metrics.yaml'))


def test_pipeline_multi_file_workflow(temp_workspace):
    """Test the pipeline command processing multiple files"""
    # Run pipeline for multiple samples
    pipeline_args = [
        'pipeline',
        os.path.join(temp_workspace, 'grnaquery.cfg'),
        '--input-dir', os.path.join(temp_workspace, 'reads'),
        '--output-dir', os.path.join(temp_workspace, 'results'),
        '--barcode-columns', 'UMI_5prime,UMI_3prime',
        '--key-column', 'protospacer',
        '--mismatches', '2',
        '--method', 'directional',
        '--metrics'
    ]
    cli = Cli(pipeline_args)
    cli.run()

    print(os.listdir(os.path.join(temp_workspace, 'results/extracted')))

    # Verify outputs for all samples
    expected_files = [
        '409-4_S1_L002_R1_001p.csv',
        '2-G1L9-M1_S9_L001_R1_001p.csv',
        '2-G1L9-M2_S12_L001_R1_001p.csv'
    ]
    
    for filename in expected_files:
        assert os.path.exists(os.path.join(temp_workspace, f'results/extracted/{filename}'))
        assert os.path.exists(os.path.join(temp_workspace, f'results/collapsed/{filename}'))
        assert os.path.exists(os.path.join(temp_workspace, f'results/counted/{filename}'))
    
    # Verify metrics files
    assert os.path.exists(os.path.join(temp_workspace, 'results/collapsed/collapse_metrics.yaml'))
    assert os.path.exists(os.path.join(temp_workspace, 'results/counted/count_metrics.yaml'))

def test_full_workflow_with_config(temp_workspace):
    """Test the full workflow using a master TOML config file"""
    # Create master config file
    config_file = os.path.join(temp_workspace, 'master_config.toml')
    with open(config_file, 'w') as f:
        f.write("""[findseq]
# No specific config needed for findseq

[collapse]
mismatches = 2
method = "directional"
sep = ","
columns = "UMI_5prime,UMI_3prime"

[count]
barcode_column = "UMI_5prime_UMI_3prime_corrected"
key_column = "protospacer"
sep = ","
detailed = false
""")

    # Step 1: Run findseq commands for each sample
    pairs = [('reads/409-4_S1_L002_R1_001.fastq.gz', 'reads/409-4_S1_L002_R2_001.fastq.gz', 'shuffle'),
             ('reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz', 'reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz', 'M1-lib'),
             ('reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz', 'reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz', 'M2-lib')]

    for read1, read2, output_name in pairs:
        findseq_args = [
            'findseq',
            os.path.join(temp_workspace, 'grnaquery.cfg'),
            '-1', os.path.join(temp_workspace, read1),
            '-2', os.path.join(temp_workspace, read2),
            '-o', os.path.join(temp_workspace, f'results/extracted/{output_name}.csv')
        ]
        cli = Cli(findseq_args)
        cli.run()

        # Verify findseq output
        assert os.path.exists(os.path.join(temp_workspace, f'results/extracted/{output_name}.csv'))

    # Step 2: Run collapse command with config
    collapse_args = [
        'collapse',
        '--config', config_file,
        '--input-dir', os.path.join(temp_workspace, 'results/extracted'),
        '--output-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--columns', 'UMI_5prime,UMI_3prime'  # Still required as it's a required argument
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'results/collapsed/{output_name}.csv'))

    # Step 3: Run count command with config
    count_args = [
        'count',
        '--config', config_file,
        '--input-dir', os.path.join(temp_workspace, 'results/collapsed'),
        '--output-dir', os.path.join(temp_workspace, 'results/counted'),
        '--barcode-column', 'UMI_5prime_UMI_3prime_corrected',  # Still required as it's a required argument
        '--key-column', 'protospacer'  # Still required as it's a required argument
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ['shuffle', 'M1-lib', 'M2-lib']:
        assert os.path.exists(os.path.join(temp_workspace, f'results/counted/{output_name}.csv'))

    # Verify config values were applied correctly
    # For collapse
    collapse_cli = Cli(collapse_args)
    assert collapse_cli.args.mismatches == 2
    assert collapse_cli.args.method == 'directional'
    assert collapse_cli.args.sep == ','

    # For count
    count_cli = Cli(count_args)
    assert count_cli.args.sep == ','
    assert count_cli.args.detailed is False

def test_pipeline_multi_file_workflow_with_config(temp_workspace):
    """Test the pipeline command processing multiple files using a master TOML config file"""
    # Create master config file
    config_file = os.path.join(temp_workspace, 'master_config.toml')
    with open(config_file, 'w') as f:
        f.write("""[pipeline]
mismatches = 2
method = "directional"
barcode_columns = "UMI_5prime,UMI_3prime"
key_column = "protospacer"
sep = ","
metrics = true

[collapse]
mismatches = 2
method = "directional"
sep = ","
columns = "UMI_5prime,UMI_3prime"

[count]
barcode_column = "UMI_5prime_UMI_3prime_corrected"
key_column = "protospacer"
sep = ","
detailed = false
""")

    # Run pipeline for multiple samples
    pipeline_args = [
        'pipeline',
        os.path.join(temp_workspace, 'grnaquery.cfg'),
        '--config', config_file,
        '--input-dir', os.path.join(temp_workspace, 'reads'),
        '--output-dir', os.path.join(temp_workspace, 'results')
    ]
    cli = Cli(pipeline_args)
    cli.run()

    # Verify config values were applied correctly
    assert cli.args.mismatches == 2
    assert cli.args.method == 'directional'
    assert cli.args.barcode_columns == 'UMI_5prime,UMI_3prime'
    assert cli.args.key_column == 'protospacer'
    assert cli.args.sep == ','
    assert cli.args.metrics is True

    # Verify outputs for all samples
    expected_files = [
        '409-4_S1_L002_R1_001p.csv',
        '2-G1L9-M1_S9_L001_R1_001p.csv',
        '2-G1L9-M2_S12_L001_R1_001p.csv'
    ]
    
    for filename in expected_files:
        assert os.path.exists(os.path.join(temp_workspace, f'results/extracted/{filename}'))
        assert os.path.exists(os.path.join(temp_workspace, f'results/collapsed/{filename}'))
        assert os.path.exists(os.path.join(temp_workspace, f'results/counted/{filename}'))
    
    # Verify metrics files
    assert os.path.exists(os.path.join(temp_workspace, 'results/collapsed/collapse_metrics.yaml'))
    assert os.path.exists(os.path.join(temp_workspace, 'results/counted/count_metrics.yaml'))

