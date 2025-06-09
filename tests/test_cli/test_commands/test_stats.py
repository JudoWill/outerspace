"""Tests for the stats command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli

def test_stats_initialization():
    """Test that stats command initializes correctly"""
    args = [
        'stats',
        'test_input.csv',
        '--umi-column', 'counts'
    ]
    cli = Cli(args)
    assert cli.args.input_files == ['test_input.csv']
    assert cli.args.umi_column == 'counts'
    assert cli.args.sep == ','
    assert cli.args.scale is None
    assert cli.args.count_column is None
    assert cli.args.allowed_list is None

def test_stats_missing_input_files():
    """Test that stats command handles missing input files"""
    args = [
        'stats',
        '--umi-column', 'counts'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_stats_missing_umi_column():
    """Test that stats command handles missing UMI column"""
    args = [
        'stats',
        'test_input.csv'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_stats_invalid_scale():
    """Test that stats command handles invalid scale value"""
    args = [
        'stats',
        'test_input.csv',
        '--umi-column', 'counts',
        '--scale', '-1'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_stats_nonexistent_input_file():
    """Test that stats command handles nonexistent input file"""
    args = [
        'stats',
        'nonexistent.csv',
        '--umi-column', 'counts'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_stats_nonexistent_umi_column():
    """Test that stats command handles nonexistent UMI column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'stats',
            temp_file.name,
            '--umi-column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_stats_nonexistent_count_column():
    """Test that stats command handles nonexistent count column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'stats',
            temp_file.name,
            '--umi-column', 'header1',
            '--count-column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_stats_nonexistent_allowed_list():
    """Test that stats command handles nonexistent allowed list file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'stats',
            temp_file.name,
            '--umi-column', 'header1',
            '--allowed-list', 'nonexistent.txt'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_stats_with_config_file():
    """Test that stats command loads and uses config file correctly"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[stats]
sep = ";"
scale = 2.0
""")
        
        # Create test input file
        input_file = os.path.join(temp_dir, 'test.csv')
        with open(input_file, 'w') as f:
            f.write('value;count\n')
            f.write('A;1\n')
            f.write('B;2\n')
            f.write('C;3\n')
        
        args = [
            'stats',
            '--config', config_file,
            input_file,
            '--umi-column', 'value',
            '--count-column', 'count'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify config values were applied
        assert cli.args.sep == ';'
        assert cli.args.scale == 2.0

def test_stats_config_override():
    """Test that command line arguments override config file settings"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[stats]
sep = ";"
scale = 2.0
""")
        
        # Create test input file
        input_file = os.path.join(temp_dir, 'test.csv')
        with open(input_file, 'w') as f:
            f.write('value|count\n')
            f.write('A|1\n')
            f.write('B|2\n')
            f.write('C|3\n')
        
        args = [
            'stats',
            '--config', config_file,
            input_file,
            '--umi-column', 'value',
            '--count-column', 'count',
            '--sep', '|',  # Override config
            '--scale', '1.5'  # Override config
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify command line args took precedence
        assert cli.args.sep == '|'
        assert cli.args.scale == 1.5

def test_stats_nonexistent_config():
    """Test that stats command handles nonexistent config file"""
    args = [
        'stats',
        '--config', 'nonexistent.toml',
        'test_input.csv',
        '--umi-column', 'value'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Configuration file not found"):
        cli.run()

def test_stats_invalid_config_section():
    """Test that stats command handles config file with invalid section"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with wrong section name
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[wrong_section]
sep = ";"
""")
        
        args = [
            'stats',
            '--config', config_file,
            'test_input.csv',
            '--umi-column', 'value'
        ]
        cli = Cli(args)
        # Should still work, just using defaults for stats-specific settings
        assert cli.args.sep == ','  # Default value

def test_stats_multiple_input_files():
    """Test that stats command handles multiple input files"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input files
        input_files = []
        for i in range(3):
            input_file = os.path.join(temp_dir, f'test{i}.csv')
            with open(input_file, 'w') as f:
                f.write('value,count\n')
                f.write('A,1\n')
                f.write('B,2\n')
                f.write('C,3\n')
            input_files.append(input_file)
        
        args = [
            'stats',
            '--umi-column', 'value',
            '--count-column', 'count'
        ] + input_files
        
        cli = Cli(args)
        cli.run()
        
        # Verify all files were processed
        assert len(cli.args.input_files) == 3 