"""Tests for the collapse command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli

def test_collapse_initialization():
    """Test that collapse command initializes correctly"""
    args = [
        'collapse',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(args)
    assert cli.args.input_dir == 'test_input'
    assert cli.args.output_dir == 'test_output'
    assert cli.args.columns == 'umi3,umi5'
    assert cli.args.mismatches == 2
    assert cli.args.method == 'directional'
    assert cli.args.sep == ','  # default value

def test_collapse_single_file_initialization():
    """Test that collapse command initializes correctly in single file mode"""
    args = [
        'collapse',
        '--input-file', 'test_input.csv',
        '--output-file', 'test_output.csv',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(args)
    assert cli.args.input_file == 'test_input.csv'
    assert cli.args.output_file == 'test_output.csv'
    assert cli.args.columns == 'umi3,umi5'
    assert cli.args.mismatches == 2
    assert cli.args.method == 'directional'
    assert cli.args.sep == ','  # default value

def test_collapse_missing_input():
    """Test that collapse command handles missing input"""
    args = [
        'collapse',
        '--output-dir', 'test_output',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_collapse_missing_output():
    """Test that collapse command handles missing output"""
    args = [
        'collapse',
        '--input-dir', 'test_input',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_collapse_missing_columns():
    """Test that collapse command handles missing columns"""
    args = [
        'collapse',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_collapse_invalid_method():
    """Test that collapse command handles invalid clustering method"""
    args = [
        'collapse',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'invalid_method'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_collapse_nonexistent_input_dir():
    """Test that collapse command handles nonexistent input directory"""
    args = [
        'collapse',
        '--input-dir', 'nonexistent_dir',
        '--output-dir', 'test_output',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_collapse_empty_input_dir():
    """Test that collapse command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            'collapse',
            '--input-dir', temp_dir,
            '--output-dir', 'test_output',
            '--columns', 'umi3,umi5',
            '--mismatches', '2',
            '--method', 'directional'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_collapse_parse_columns():
    """Test that column parsing works correctly"""
    args = [
        'collapse',
        '--input-dir', 'test_input',
        '--output-dir', 'test_output',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(args)
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ['umi3', 'umi5']
    
    # Test with spaces
    cli.args.columns = 'umi3, umi5'
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ['umi3', 'umi5']
    
    # Test with single column
    cli.args.columns = 'umi3'
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ['umi3']

def test_collapse_single_file_missing_output_file():
    """Test that collapse command handles missing output file for single file mode"""
    args = [
        'collapse',
        '--input-file', 'test_input.csv',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_collapse_single_file_nonexistent_input():
    """Test that collapse command handles nonexistent input file in single file mode"""
    args = [
        'collapse',
        '--input-file', 'nonexistent.csv',
        '--output-file', 'test_output.csv',
        '--columns', 'umi3,umi5',
        '--mismatches', '2',
        '--method', 'directional'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_collapse_single_file_basic():
    """Test basic single file collapse functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, 'test.csv')
        with open(input_file, 'w') as f:
            f.write('umi3,umi5\n')
            f.write('AAA,CCC\n')
            f.write('AAT,CCC\n')  # One mismatch from AAA
            f.write('GGG,TTT\n')
        
        # Create output file path
        output_file = os.path.join(temp_dir, 'output.csv')
        
        args = [
            'collapse',
            '--input-file', input_file,
            '--output-file', output_file,
            '--columns', 'umi3,umi5',
            '--mismatches', '1',
            '--method', 'directional'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, 'r') as f:
            content = f.read()
            assert 'umi3,umi5,umi3_umi5_corrected' in content
            # The first two rows should have the same corrected barcode
            assert 'AAA,CCC' in content
            assert 'AAT,CCC' in content
            assert 'GGG,TTT' in content

