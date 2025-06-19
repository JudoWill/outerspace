"""Tests for the gini command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli

def test_gini_initialization():
    """Test that gini command initializes correctly"""
    args = [
        'gini',
        'test_input.csv',
        '--column', 'counts'
    ]
    cli = Cli(args)
    assert cli.args.input_file == 'test_input.csv'
    assert cli.args.column == 'counts'
    assert cli.args.sep == ','
    assert cli.args.scale is None
    assert cli.args.count_column is None
    assert cli.args.allowed_list is None

def test_gini_missing_input_file():
    """Test that gini command handles missing input file"""
    args = [
        'gini',
        '--column', 'counts'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_gini_missing_column():
    """Test that gini command handles missing column"""
    args = [
        'gini',
        'test_input.csv'
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2

def test_gini_invalid_scale():
    """Test that gini command handles invalid scale value"""
    args = [
        'gini',
        'test_input.csv',
        '--column', 'counts',
        '--scale', '-1'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_gini_nonexistent_input_file():
    """Test that gini command handles nonexistent input file"""
    args = [
        'gini',
        'nonexistent.csv',
        '--column', 'counts'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()

def test_gini_nonexistent_column():
    """Test that gini command handles nonexistent column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_gini_nonexistent_count_column():
    """Test that gini command handles nonexistent count column"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'header1',
            '--count-column', 'nonexistent'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_gini_nonexistent_allowed_list():
    """Test that gini command handles nonexistent allowed list file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv') as temp_file:
        temp_file.write('header1,header2\nvalue1,value2\n')
        temp_file.flush()
        
        args = [
            'gini',
            temp_file.name,
            '--column', 'header1',
            '--allowed-list', 'nonexistent.txt'
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()

def test_gini_with_config_file():
    """Test that gini command loads and uses config file correctly"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[gini]
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
            'gini',
            '--config', config_file,
            input_file,
            '--column', 'value',
            '--count-column', 'count'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify config values were applied
        assert cli.args.sep == ';'
        assert cli.args.scale == 2.0

def test_gini_config_override():
    """Test that command line arguments override config file settings"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[gini]
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
            'gini',
            '--config', config_file,
            input_file,
            '--column', 'value',
            '--count-column', 'count',
            '--sep', '|',  # Override config
            '--scale', '1.5'  # Override config
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify command line args took precedence
        assert cli.args.sep == '|'
        assert cli.args.scale == 1.5

def test_gini_nonexistent_config():
    """Test that gini command handles nonexistent config file"""
    args = [
        'gini',
        '--config', 'nonexistent.toml',
        'test_input.csv',
        '--column', 'value'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Configuration file not found"):
        cli.run()

def test_gini_invalid_config_section():
    """Test that gini command handles config file with invalid section"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with wrong section name
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[wrong_section]
sep = ";"
""")
        
        args = [
            'gini',
            '--config', config_file,
            'test_input.csv',
            '--column', 'value'
        ]
        cli = Cli(args)
        # Should still work, just using defaults for gini-specific settings
        assert cli.args.sep == ','  # Default value
