"""Tests for the count command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli

def test_count_initialization():
    """Test that count command initializes correctly"""
    args = [
        'count',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(args)
    assert cli.args.input_dir == 'test_input'
    assert cli.args.output_dir == 'test_output'
    assert cli.args.barcode_column == 'umi3_umi5_corrected'
    assert cli.args.key_column == 'protospacer'
    assert cli.args.sep == ','
    assert cli.args.detailed is False

def test_count_single_file_initialization():
    """Test that count command initializes correctly in single file mode"""
    args = [
        'count',
        '--input-file', 'test_input.csv',
        '--output-file', 'test_output.csv',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(args)
    assert cli.args.input_file == 'test_input.csv'
    assert cli.args.output_file == 'test_output.csv'
    assert cli.args.barcode_column == 'umi3_umi5_corrected'
    assert cli.args.key_column == 'protospacer'
    assert cli.args.sep == ','
    assert cli.args.detailed is False

def test_count_missing_input():
    """Test that count command handles missing input"""
    args = [
        'count',
        '--output-dir', 'test_output',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_count_missing_output():
    """Test that count command handles missing output"""
    args = [
        'count',
        '--input-dir', 'test_input',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_count_missing_barcode_column():
    """Test that count command handles missing barcode column"""
    args = [
        'count',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--key-column', 'protospacer'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_count_missing_key_column():
    """Test that count command handles missing key column"""
    args = [
        'count',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--barcode-column', 'umi3_umi5_corrected'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_count_invalid_downsample():
    """Test that count command handles invalid downsample value"""
    args = [
        'count',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer',
        '--downsample', '-1'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_count_nonexistent_input_dir():
    """Test that count command handles nonexistent input directory"""
    args = [
        'count',
        '--input-dir', 'nonexistent_dir',
        '--output-dir', 'test_output',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_count_empty_input_dir():
    """Test that count command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            'count',
            '--input-dir', temp_dir,
            '--output-dir', 'test_output',
            '--barcode-column', 'umi3_umi5_corrected',
            '--key-column', 'protospacer'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_count_single_file_missing_output_file():
    """Test that count command handles missing output file for single file mode"""
    args = [
        'count',
        '--input-file', 'test_input.csv',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_count_single_file_nonexistent_input():
    """Test that count command handles nonexistent input file in single file mode"""
    args = [
        'count',
        '--input-file', 'nonexistent.csv',
        '--output-file', 'test_output.csv',
        '--barcode-column', 'umi3_umi5_corrected',
        '--key-column', 'protospacer'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_count_single_file_basic():
    """Test basic single file counting functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, 'test.csv')
        with open(input_file, 'w') as f:
            f.write('protospacer,umi3_umi5_corrected\n')
            f.write('key1,bar1\n')
            f.write('key1,bar2\n')
            f.write('key2,bar3\n')
        
        # Create output file path
        output_file = os.path.join(temp_dir, 'output.csv')
        
        args = [
            'count',
            '--input-file', input_file,
            '--output-file', output_file,
            '--barcode-column', 'umi3_umi5_corrected',
            '--key-column', 'protospacer'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, 'r') as f:
            content = f.read()
            assert 'key1,2' in content
            assert 'key2,1' in content

def test_count_single_file_with_allowed_list():
    """Test single file counting with allowed list"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, 'test.csv')
        with open(input_file, 'w') as f:
            f.write('protospacer,umi3_umi5_corrected\n')
            f.write('key1,bar1\n')
            f.write('key1,bar2\n')
            f.write('key2,bar3\n')
            f.write('key3,bar4\n')  # This key should be filtered out
        
        # Create allowed list file
        allowed_list = os.path.join(temp_dir, 'allowed.txt')
        with open(allowed_list, 'w') as f:
            f.write('key1\n')
            f.write('key2\n')
        
        # Create output file path
        output_file = os.path.join(temp_dir, 'output.csv')
        
        args = [
            'count',
            '--input-file', input_file,
            '--output-file', output_file,
            '--barcode-column', 'umi3_umi5_corrected',
            '--key-column', 'protospacer',
            '--allowed-list', allowed_list
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, 'r') as f:
            content = f.read()
            assert 'key1,2' in content
            assert 'key2,1' in content
            assert 'key3' not in content

